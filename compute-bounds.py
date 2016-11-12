#!/usr/bin/python

import sys
import argparse

try:
    from gurobipy import *
except ImportError:
    raise Exception('gurobipy not found.\nthis code requires gurobi to be installed in python')

from parse_matpower import *
from qc_lib import *
import time

version = '1.0.0'

def build_parser():
    parser = argparse.ArgumentParser(
        description='''compute-bounds is a python/gurobi based script for 
            computing tight bounds on voltage magnitudes and line phase angle 
            differences in matpower Optimal Power Flow datasets.  
            The theory behind this code can be found in 'Strengthening 
            Convex Relaxations with Bound Tightening for Power Network 
            Optimization', Carleton Coffrin, Hassan L. Hijazi, and 
            Pascal Van Hentenryck.''',

        epilog='''Please file bugs at https://github.com/ccoffrin/ac-opf-bounds''',
    )
    parser.add_argument('file', help='a matpower case file to process (.m)')
    parser.add_argument('-loa', '--linear_oa', help='use a linear outer-approximation of the QC model (more scalable)', action='store_true')
    parser.add_argument('-i', '--iterations', help='limit the number of iterations (0 computes the base relaxation without tightening)', type=int, default=sys.maxint)
    parser.add_argument('-l', '--large', help='configures the algorithm for processing cases with over 1000 buses', action='store_true')
    parser.add_argument('-pad', help='override test case phase angle difference bounds with given value (in radians)', type=float)

    # rough output levels
    # 0 minimal
    # 1-3 tbd
    # 4 add all algorithm output
    # 5 add qc-model output
    # 6 add gurobi solve output
    parser.add_argument('-o', '--output', help='controls the output level (0-6)', type=int, default=0)
    parser.add_argument('-v', '--version', action='version', version=version)

    return parser


def main(args):
    '''reads a matpower case files and processes them based on command 
    line arguments.

    Args:
        args: an argparse data structure
    '''

    case = parse_mp_case(args.file)
    case = case.make_per_unit()
    case = case.make_radians()
    if args.pad != None:
        if args.pad <= 0:
            raise Exception('pad parameter should be a positive number')
        case = case.replace_phase_angle_difference(args.pad)
    case = case.update_line_limits()

    ### WARNING ###
    # for small cases (less than 1000 buses) tighter values can be used
    # these parameters have been tuned to be robust on cases with over 1000 buses
    #

    # optimality tolerance of the solver
    opt_tol = 1e-6

    # the number significant figures when tightening bounds from and objective value
    # note, should be less than opt_tol
    output_tol = 1e+5

    # stop tightening a variable's bounds when the range is smaller than this value
    min_bound_width  = 1e-3

    # tolerance for determining when a fixpoint is reached
    improvement_tolerace = 1e-3

    time_0 = time.time()

    m = Model('QC')
    if args.output <= 5:
        m.setParam('OutputFlag', 0)
    m.setParam('LogFile', '')
    m.setParam('Crossover', 0)
    m.setParam('OptimalityTol', opt_tol)
    m.setParam('FeasibilityTol', 1e-8)

    # special parameters that improve reliability on large cases
    if args.large:
        m.setParam('Method', 2)
        m.setParam('TimeLimit', 20) # seconds
        m.setParam('BarIterLimit', 200)


    qc_model = QCModel(m, case, output_level=args.output)

    if args.linear_oa or args.large:
        qc_model.build_linear()
    else:
        qc_model.build()

    print('')
    print('Initial Relaxation:')
    print(m)

    obj = qc_model.setMinCostObjective()
    m.optimize()

    print('Status: %s' % m.status)

    obj_val = float('inf')
    if m.status == 2:
        obj_val = obj.getValue()
        print('Objective: %f' % obj_val)
    else:
        raise Exception('relaxation at the root node failed... \nmost likely the case is infeasible.')


    base_obj_val = obj_val

    #print('INFO_1, %s, %s, %s, %s' % (args.file, m.status, m.runtime, obj_val))

    v = qc_model.v
    td = qc_model.td

    b_count = len(case.bus)
    l_count = len(case.branch)
    bp_count = len(qc_model.bus_pairs)
    iteration = 0
    iter_time = time.time()
    max_subproblem_runtime = 0
    total_parallel_runtime = m.runtime
    total_v_reduction = b_count
    total_td_reduction = bp_count
    max_v_reduction = 1
    max_td_reduction = 1

    v_lb_init = {i:x.lb for i,x in v.iteritems()}
    v_ub_init = {i:x.ub for i,x in v.iteritems()}
    td_lb_init = {i:x.lb for i,x in td.iteritems()}
    td_ub_init = {i:x.ub for i,x in td.iteritems()}

    print('')
    print('id, name, buses, lines, obj, model time(sec.), iter., time, max prob. time, sum v rng, sum pad rng, avg. v rng, avg. pad rng, avg. v reduc., avg. pad reduc., pad sign fixed, model status')
    #TODO turn up this tolerance to 0.0001, will improve some nesta_case30_ieee.m at least
    while (total_td_reduction/bp_count > improvement_tolerace or total_v_reduction/b_count > improvement_tolerace) and args.iterations > iteration:
    #while (max_td_reduction > improvement_tolerace or max_td_reduction > improvement_tolerace) and args.iterations > iteration:

        v_range = sum([x.ub-x.lb for x in v.values()])
        td_range = sum([x.ub-x.lb for x in td.values()])
        dir_count = sum([x.ub <= 0 or x.lb >= 0 for x in td.values()])

        print('ITER_DATA, %s, %d, %d, %f, %f, %d, %f, %f, %f, %f, %f, %f, %f, %f, %d, %d' % (case.name, b_count, l_count, obj_val, m.runtime, iteration, time.time() - iter_time, max_subproblem_runtime, v_range, td_range, v_range/b_count, td_range/bp_count, total_v_reduction/b_count, total_td_reduction/bp_count, dir_count, m.status))

        max_subproblem_runtime = 0
        total_v_reduction = 0
        total_td_reduction = 0
        max_v_reduction = 0
        max_td_reduction = 0
        iter_time = time.time()

        if args.output >= 4:
            sys.stdout.write('v ')
        for i,x in v.iteritems():
            lb = None
            qc_model.setObjective(x,GRB.MINIMIZE)
            m.optimize()
            if args.output >= 4:
                sys.stdout.write(' '+str(m.status))
                sys.stdout.flush()
            if m.status == 2:
                nlb = floor(output_tol*x.x)/output_tol
                if nlb > x.lb:
                    lb = nlb

            max_subproblem_runtime = max(max_subproblem_runtime, m.runtime)

            ub = None
            qc_model.setObjective(x,GRB.MAXIMIZE)
            m.optimize()
            if args.output >= 4:
                sys.stdout.write(' '+str(m.status))
                sys.stdout.flush()
            if m.status == 2:
                nub = ceil(output_tol*x.x)/output_tol
                if nub < x.ub:
                    ub = nub

            max_subproblem_runtime = max(max_subproblem_runtime, m.runtime)

            #print('%s, %f, %f, %f, %f, %f' % (i, x.lb, x.ub, lb, ub, qc_model.m.runtime))
            if lb != None and lb > x.ub:
                print('\nwarning: on %s computed lb %f was above assigned ub %f' % (x, lb, x.ub))
                lb = x.lb
            if ub != None and ub < x.lb:
                print('\nwarning: on %s computed ub %f is below assigned lb %f' % (x, ub, x.lb))
                ub = x.ub
                
            if lb == None:
                lb = x.lb
            if ub == None:
                ub = x.ub

            if ub - lb >= min_bound_width:
                v_reduction = (x.ub - x.lb) - (ub - lb)
                x.ub = ub
                x.lb = lb 
                qc_model.boundChanged(x)
            elif x.ub - lb >= min_bound_width:
                v_reduction = (x.ub - x.lb) - (x.ub - lb)
                x.lb = lb 
                qc_model.boundChanged(x)
            elif ub - x.lb >= min_bound_width:
                v_reduction = (x.ub - x.lb) - (ub - x.lb)
                x.ub = ub 
                qc_model.boundChanged(x)
            else:
                v_reduction = 0

            total_v_reduction += v_reduction
            max_v_reduction = max(max_v_reduction, v_reduction)
        if args.output >= 4:
            sys.stdout.write('\n')

        if args.output >= 4:
            sys.stdout.write('pad ')
        for i,x in td.iteritems():
            lb = None
            qc_model.setObjective(x,GRB.MINIMIZE)
            m.optimize()
            if args.output >= 4:
                sys.stdout.write(' '+str(m.status))
                sys.stdout.flush()
            if m.status == 2:
                nlb = floor(output_tol*x.x)/output_tol
                if nlb > x.lb:
                    lb = nlb
            
            max_subproblem_runtime = max(max_subproblem_runtime, m.runtime)

            ub = None
            qc_model.setObjective(x,GRB.MAXIMIZE)
            m.optimize()
            if args.output >= 4:
                sys.stdout.write(' '+str(m.status))
                sys.stdout.flush()
            if m.status == 2:
                nub = ceil(output_tol*x.x)/output_tol
                if nub < x.ub:
                    ub = nub

            max_subproblem_runtime = max(max_subproblem_runtime, m.runtime)


            #print('%s, %f, %f, %f, %f, %f' % (i, x.lb, x.ub, lb, ub, qc_model.m.runtime))
            if lb != None and lb > x.ub:
                print('\nwarning: on %s computed lb %f was above assigned ub %f' % (x, lb, x.ub))
                lb = x.lb
            if ub != None and ub < x.lb:
                print('\nwarning: on %s computed ub %f is below assigned lb %f' % (x, ub, x.lb))
                ub = x.ub
                
            if lb == None:
                lb = x.lb
            if ub == None:
                ub = x.ub

            if ub - lb >= min_bound_width:
                td_reduction = (x.ub - x.lb) - (ub - lb)
                x.ub = ub
                x.lb = lb 
                qc_model.boundChanged(x)
            elif x.ub - lb >= min_bound_width:
                td_reduction = (x.ub - x.lb) - (x.ub - lb)
                x.lb = lb 
                qc_model.boundChanged(x)
            elif ub - x.lb >= min_bound_width:
                td_reduction = (x.ub - x.lb) - (ub - x.lb)
                x.ub = ub 
                qc_model.boundChanged(x)
            else:
                td_reduction = 0

            total_td_reduction += td_reduction
            max_td_reduction = max(max_td_reduction, td_reduction)
        if args.output >= 4:
            sys.stdout.write('\n')

        update_time = time.time()
        qc_model.update()
        #qc_model.update_all()

        #v_bound_ranges = [x.ub - x.lb for i,x in v.iteritems()]
        #td_bound_ranges = [x.ub - x.lb for i,x in td.iteritems()]
        #print('v range: ', min(v_bound_ranges), max(v_bound_ranges))
        #print('td range:', min(td_bound_ranges), max(td_bound_ranges))


        obj = qc_model.setMinCostObjective()
        #m.setParam('OutputFlag',1)
        m.optimize()
        #m.setParam('OutputFlag',0)
        #print('solve status: ', m.status)
        if m.status == 2:
            obj_val = obj.getValue()
        #if m.status == 13:
        #    m.setParam('NumericFocus',3)

        total_parallel_runtime += max_subproblem_runtime
        iteration += 1

    v_range = sum([x.ub-x.lb for x in v.values()])
    td_range = sum([x.ub-x.lb for x in td.values()])
    dir_count = sum([x.ub <= 0 or x.lb >= 0 for x in td.values()])

    print('ITER_DATA, %s, %d, %d, %f, %f, %d, %f, %f, %f, %f, %f, %f, %f, %f, %d, %d' % (case.name, b_count, l_count, obj_val, m.runtime, iteration, time.time() - iter_time, max_subproblem_runtime, v_range, td_range, v_range/b_count, td_range/bp_count, total_v_reduction/b_count, total_td_reduction/bp_count, dir_count, m.status))

    print('')
    print('name, bus index, v min, v max')
    print('//START_BUS_DATA//')
    for i,b in qc_model.buses.iteritems():
        print('%s, %d, %f, %f' % (case.name, i, max(v_lb_init[i], floor(output_tol*v[i].lb)/output_tol), min(v_ub_init[i], ceil(output_tol*v[i].ub)/output_tol)))
    print('//END_BUS_DATA//\n')

    print('name, line index, pad min, pad max')
    print('//START_LINE_DATA//')
    for arc,br in qc_model.arcs_from.iteritems():
        bp = (arc[1],arc[2])
        print('%s, %d, %d, %d, %f, %f' % (case.name, arc[0]+1, arc[1], arc[2], max(td_lb_init[bp], floor(output_tol*td[bp].lb)/output_tol), min(td_ub_init[bp], ceil(output_tol*td[bp].ub)/output_tol)))
    print('//END_LINE_DATA//\n')

    v_range_init = sum([v_ub_init[i]-v_lb_init[i] for i in v.keys()])
    td_range_init = sum([td_ub_init[i]-td_lb_init[i] for i in td.keys()])


    print('id, name, buses, lines, start obj, end obj, time(sec.), iter., par time(sec.), avg. v rng start, avg. v rng end, avg. pad rng start, avg. pad rng end, pad sign fixed')
    print('SMRY_DATA, %s, %d, %d, %f, %f, %f, %d, %f, %f, %f, %f, %f, %d' % (case.name, b_count, l_count, base_obj_val, obj_val, time.time() - time_0, iteration, total_parallel_runtime, v_range_init/b_count, v_range/b_count, td_range_init/bp_count, td_range/bp_count, dir_count))

    #print('ALL_DATA, %s, %d, %d, %f, %d, %f, -inf' % (case.name, b_count, l_count, obj_val, 100, time.time() - time_0))

    #print('INFO_2, %s, %s, %s, %s' % (args.file, m.status, m.runtime, obj_val))

    if  m.status == 3 or m.status == 4:
        print('infeasible relaxation!!!')
        print('either there is a bug in this algorithm or there are numerical inconsistencies in the solver.')
        # m.computeIIS()
        # print('\nThe following constraint(s) cannot be satisfied:')
        # for v in m.getVars():
        #     if v.IISLB:
        #         print('ISS LB %s' % v.varName)
        #     if v.IISUB:
        #         print('ISS UB %s' % v.varName)
        #     #print v.varName, v.x
        # for c in m.getConstrs():
        #     if c.IISConstr:
        #         print('%s' % c.constrName)
        quit()
        
    if  m.status == 12:
        print('Numerical trouble in the solver. =(')
        quit()


if __name__ == '__main__':
    parser = build_parser()
    main(parser.parse_args())

