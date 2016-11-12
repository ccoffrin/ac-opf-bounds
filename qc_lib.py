from math import cos, sin, tan, pi, sqrt, ceil, floor
from gurobipy import *
import time


class CurrentSqrRelaxScheme(object):
    def __init__(self, model, x, y, s, tr, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.x = x
        self.y = y

        self.s = s
        self.tr = tr

        self.update()

    def update(self):
        self.y.lb = 0 
        self.y.ub = (self.s*self.tr/self.x.lb)**2


class CurrentSqrBoundsRelaxScheme(object):
    def __init__(self, model, wfb, wtb, wr, wi, qfb, br, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.wfb = wfb
        self.wtb = wtb
        self.wr = wr
        self.wi = wi
        self.qfb = qfb

        gl =  br.r / (br.r**2 + br.x**2)
        bl = -br.x / (br.r**2 + br.x**2)

        tr = br.ratio
        if tr == 0:
            tr = 1.0

        self.ysqr = gl**2 + bl**2
        self.tr = tr
        self.tr2 = self.tr**2
        self.cc = tr*cos(br.angle)
        self.dd = tr*sin(br.angle)
        self.ch = br.b
        self.s = br.rateA

        self.c_1 = None
        self.c_2 = None

        self.update()

    def update(self):
        if self.c_1 == None:
            self.c_1 = self.model.addConstr(self.ysqr*(self.wfb/self.tr2 + self.wtb - 2*(self.cc*self.wr + self.dd*self.wi)/self.tr2) - self.ch*self.qfb - ((self.ch/2)/self.tr)**2*self.wfb >= 0, self.prefix+'c_lb_'+self.postfix)

        if self.c_2 != None:
            self.model.remove(self.c_2)

        self.c_2 = self.model.addConstr(self.ysqr*(self.wfb/self.tr2 + self.wtb - 2*(self.cc*self.wr + self.dd*self.wi)/self.tr2) - self.ch*self.qfb - ((self.ch/2)/self.tr)**2*self.wfb <= (self.s*self.tr)**2/self.wfb.lb, self.prefix+'c_ub_'+self.postfix)


class WPADRelaxScheme(object):
    def __init__(self, model, x, y, td, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.x = x
        self.y = y
        self.td = td

        self.c_1 = None
        self.c_2 = None

        self.update()

    def update(self):
        if self.c_1 != None:
            self.model.remove(self.c_1)
        if self.c_2 != None:
            self.model.remove(self.c_2)

        self.c_1 = self.model.addConstr(self.y <= tan(self.td.ub)*self.x, self.prefix+'w_pad_ub'+self.postfix)
        self.c_2 = self.model.addConstr(self.y >= tan(self.td.lb)*self.x, self.prefix+'w_pad_lb'+self.postfix)


class SquareRelaxScheme(object):
    def __init__(self, model, x, y, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.x = x
        self.y = y

        self.c_1 = None
        self.c_2 = None

        self.update()

    def update(self):
        if self.c_1 == None:
            self.c_1 = self.model.addConstr(self.y >= self.x*self.x, self.prefix+'sqr_lb'+self.postfix)
        if self.c_2 != None:
            self.model.remove(self.c_2)

        #NOTE, turning this off seems to help reliblity on nesta_case189_edin__api.m
        #if self.x.ub - self.x.lb <= epsilon:
        #    self.c_1 = self.model.addConstr(self.y == self.x.lb**2, self.prefix+'sqr_eq'+self.postfix)
        #    return

        self.y.ub = self.x.ub**2
        self.y.lb = self.x.lb**2
        self.c_2 = self.model.addConstr(self.y <= (self.x.ub+self.x.lb)*self.x - self.x.ub*self.x.lb, self.prefix+'sqr_ub'+self.postfix)


class LinSquareRelaxScheme(object):
    def __init__(self, model, x, y, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.x = x
        self.y = y

        self.c_1 = None
        self.c_2 = None
        self.c_3 = None
        self.c_4 = None

        self.update()

    def update(self):
        if self.c_1 != None:
            self.model.remove(self.c_1)
        if self.c_2 != None:
            self.model.remove(self.c_2)
        if self.c_3 != None:
            self.model.remove(self.c_3)
        if self.c_4 != None:
            self.model.remove(self.c_4)

        self.y.ub = self.x.ub**2
        self.y.lb = self.x.lb**2
        self.c_1 = self.model.addConstr(self.y <= (self.x.ub+self.x.lb)*self.x - self.x.ub*self.x.lb, self.prefix+'sqr_ub'+self.postfix)

        #see to be causing problems
        #p1 = self.x.lb
        #self.c_2 = self.model.addConstr(self.y >= 2*p1*self.x - p1**2, self.prefix+'sqr_lb_1'+self.postfix)
        #p2 = self.x.ub
        #self.c_3 = self.model.addConstr(self.y >= 2*p2*self.x - p2**2, self.prefix+'sqr_lb_2'+self.postfix)
        
        p3 = (self.x.ub+self.x.lb)/2
        self.c_4 = self.model.addConstr(self.y >= 2*p3*self.x - p3**2, self.prefix+'sqr_lb_3'+self.postfix)


class ProdRelaxScheme(object):
    def __init__(self, model, x, y, z, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.x = x
        self.y = y
        self.z = z

        self.c_1 = None
        self.c_2 = None
        self.c_3 = None
        self.c_4 = None

        self.update()
    
    def update(self):
        if self.c_1 != None:
            self.model.remove(self.c_1)
        if self.c_2 != None:
            self.model.remove(self.c_2)
        if self.c_3 != None:
            self.model.remove(self.c_3)
        if self.c_4 != None:
            self.model.remove(self.c_4)

        xlb = self.x.lb
        xub = self.x.ub
        ylb = self.y.lb
        yub = self.y.ub

        touched = False
        if xlb >= 0 and ylb >= 0:
            self.z.lb = xlb*ylb
            self.z.ub = xub*yub
            touched = True
        if xlb >= 0 and yub <= 0:
            self.z.lb = xub*ylb
            self.z.ub = xlb*yub
            touched = True
        if xlb >= 0 and ylb < 0 and yub > 0:
            self.z.lb = xub*ylb
            self.z.ub = xub*yub
            touched = True
        assert(touched)

        self.c_1 = self.model.addConstr(self.z >= self.x.lb*self.y + self.y.lb*self.x - self.x.lb*self.y.lb, self.prefix+'mac_lb1'+self.postfix)
        self.c_2 = self.model.addConstr(self.z >= self.x.ub*self.y + self.y.ub*self.x - self.x.ub*self.y.ub, self.prefix+'mac_lb2'+self.postfix)
        self.c_3 = self.model.addConstr(self.z <= self.x.lb*self.y + self.y.ub*self.x - self.x.lb*self.y.ub, self.prefix+'mac_ub1'+self.postfix)
        self.c_4 = self.model.addConstr(self.z <= self.x.ub*self.y + self.y.lb*self.x - self.x.ub*self.y.lb, self.prefix+'mac_ub2'+self.postfix)


class SineRelaxScheme(object):
    def __init__(self, model, x, y, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.x = x
        self.y = y

        self.c_1 = None
        self.c_2 = None

        self.update()

    def update(self):
        ub = self.x.ub
        lb = self.x.lb

        if not (lb >= -pi/2 and ub <= pi/2):
            raise Exception('phase angle bounds are not within (-pi/2, pi/2) \noverride with use "-pad" argument')

        if self.c_1 != None:
            self.model.remove(self.c_1)
        if self.c_2 != None:
            self.model.remove(self.c_2)

        #if ub - lb <= epsilon:
        #    self.c_1 = self.model.addConstr(self.y == sin(lb), self.prefix+'sin_eq'+self.postfix)
        #    return

        self.y.ub = sin(ub)
        self.y.lb = sin(lb)

        max_ad = max(abs(lb),abs(ub))
        if lb < 0 and ub > 0:
            self.c_1 = self.model.addConstr(self.y <= cos(max_ad/2)*(self.x - max_ad/2) + sin(max_ad/2), self.prefix+'sin_ub'+self.postfix)
            self.c_2 = self.model.addConstr(self.y >= cos(max_ad/2)*(self.x + max_ad/2) - sin(max_ad/2), self.prefix+'sin_lb'+self.postfix)
        if ub <= 0:
            self.c_1 = self.model.addConstr(self.y <= (sin(lb) - sin(ub))/(lb-ub)*(self.x - lb) + sin(lb), self.prefix+'sin_ub'+self.postfix)
            self.c_2 = self.model.addConstr(self.y >= cos(max_ad/2)*(self.x + max_ad/2) - sin(max_ad/2), self.prefix+'sin_lb'+self.postfix)
        if lb >= 0:
            self.c_1 = self.model.addConstr(self.y <= cos(max_ad/2)*(self.x - max_ad/2) + sin(max_ad/2), self.prefix+'sin_ub'+self.postfix)
            self.c_2 = self.model.addConstr(self.y >= (sin(lb) - sin(ub))/(lb-ub)*(self.x - lb) + sin(lb), self.prefix+'sin_lb'+self.postfix)


class LinCircleRelaxScheme(object):
    def __init__(self, model, x, y, s, cuts=4, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.x = x
        self.y = y
        self.s = s

        self.cuts = range(0,cuts)

        step = 2*pi/len(self.cuts)
        offset = step/2
        self.rads = [i*step+offset for i in self.cuts] 
        self.c_list = [None for i in self.cuts]
        #print step, offset, len(self.cuts), self.rads
        
        self.update()

    def update(self):
        for c in self.cuts:
            if self.c_list[c] == None:
                px = self.s*cos(self.rads[c])
                py = self.s*sin(self.rads[c])
                #print px, py
                self.c_list[c] = self.model.addConstr(0 >= 2*px*(self.x - px) + 2*py*(self.y - py), self.prefix+'circ_ub_'+str(c+1)+self.postfix)


class CosineRelaxScheme(object):
    def __init__(self, model, x, y, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.x = x
        self.y = y

        self.c_1 = None
        self.c_2 = None

        self.update()

    def update(self):
        ub = self.x.ub
        lb = self.x.lb

        if not (lb >= -pi/2 and ub <= pi/2):
            raise Exception('phase angle bounds are not within (-pi/2, pi/2) \noverride with the "-pad" argument')

        if self.c_1 != None:
            self.model.remove(self.c_1)
        if self.c_2 != None:
            self.model.remove(self.c_2)

        if lb < 0 and ub > 0:
            self.y.ub = 1.0
            self.y.lb = min(cos(lb),cos(ub))
        if ub <= 0:
            self.y.ub = cos(ub)
            self.y.lb = cos(lb)
        if lb >= 0:
            self.y.ub = cos(lb)
            self.y.lb = cos(ub)

        max_ad = max(abs(lb),abs(ub))
        self.c_1 = self.model.addConstr(self.y <= 1 - (1-cos(max_ad))/(max_ad*max_ad)*(self.x*self.x), self.prefix+'cos_ub'+self.postfix)
        self.c_2 = self.model.addConstr(self.y >= (cos(lb) - cos(ub))/(lb-ub)*(self.x - lb) + cos(lb), self.prefix+'cos_lb'+self.postfix)


class LinCosineRelaxScheme(object):
    def __init__(self, model, x, y, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.x = x
        self.y = y

        self.c_1 = None
        self.c_2 = None
        self.c_3 = None
        self.c_4 = None

        self.update()

    def update(self):
        ub = self.x.ub
        lb = self.x.lb

        if not (lb >= -pi/2 and ub <= pi/2):
            raise Exception('phase angle bounds are not within (-pi/2, pi/2) \noverride with the "-pad" argument')

        if self.c_1 != None:
            self.model.remove(self.c_1)
        if self.c_2 != None:
            self.model.remove(self.c_2)
        if self.c_3 != None:
            self.model.remove(self.c_3)
        if self.c_4 != None:
            self.model.remove(self.c_4)

        #if ub - lb <= epsilon:
        #    self.c_1 = self.model.addConstr(self.y == cos(lb), self.prefix+'cos_eq'+self.postfix)
        #    return

        if lb < 0 and ub > 0:
            self.y.ub = 1.0
            self.y.lb = min(cos(lb),cos(ub))
        if ub <= 0:
            self.y.ub = cos(ub)
            self.y.lb = cos(lb)
        if lb >= 0:
            self.y.ub = cos(lb)
            self.y.lb = cos(ub)

        max_ad = max(abs(lb),abs(ub))
        self.c_1 = self.model.addConstr(self.y >= (cos(lb) - cos(ub))/(lb-ub)*(self.x - lb) + cos(lb), self.prefix+'cos_lb'+self.postfix)

        #p1 = self.x.lb
        #self.c_2 = self.model.addConstr(self.y <= -sin(p1)*(self.x - p1) + cos(p1), self.prefix+'cos_ub_1'+self.postfix)
        #p2 = self.x.ub
        #self.c_3 = self.model.addConstr(self.y <= -sin(p2)*(self.x - p2) + cos(p2), self.prefix+'cos_ub_2'+self.postfix)
        p3 = (self.x.ub+self.x.lb)/2
        self.c_4 = self.model.addConstr(self.y <= -sin(p3)*(self.x - p3) + cos(p3), self.prefix+'cos_ub_3'+self.postfix)


class LNCScheme(object):
    def __init__(self, model, wf, wt, wr, wi, vf, vt, td, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.wf = wf
        self.wt = wt
        self.wr = wr
        self.wi = wi

        self.vf = vf
        self.vt = vt
        self.td = td

        self.c_1 = None
        self.c_2 = None

        self.update()

    def update(self):
        vfub = self.vf.ub
        vflb = self.vf.lb
        vtub = self.vt.ub
        vtlb = self.vt.lb
        tdub = self.td.ub
        tdlb = self.td.lb

        if not (tdlb >= -pi/2 and tdub <= pi/2):
            raise Exception('phase angle bounds are not within (-pi/2, pi/2) \noverride with the "-pad" argument')

        if self.c_1 != None:
            self.model.remove(self.c_1)
        if self.c_2 != None:
            self.model.remove(self.c_2)

        phi = (tdub + tdlb)/2
        d   = (tdub - tdlb)/2

        sf = vflb + vfub
        st = vtlb + vtub

        self.c_1 = self.model.addConstr(sf*st*(cos(phi)*self.wr + sin(phi)*self.wi) - vtub*cos(d)*st*self.wf - vfub*cos(d)*sf*self.wt >=  vfub*vtub*cos(d)*(vflb*vtlb - vfub*vtub), self.prefix+'lnc_ub'+self.postfix)
        self.c_2 = self.model.addConstr(sf*st*(cos(phi)*self.wr + sin(phi)*self.wi) - vtlb*cos(d)*st*self.wf - vflb*cos(d)*sf*self.wt >= -vflb*vtlb*cos(d)*(vflb*vtlb - vfub*vtub), self.prefix+'lnc_lb'+self.postfix)


class LinSOCRelaxScheme(object):
    # relaxation for the equation a^2 + b^2 - c*d <= 0
    def __init__(self, model, a, b, c, d, prefix='', postfix=''):
        print('Warning:This constraints is incomplete!')
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.a = a
        self.b = b
        self.c = c
        self.d = d

        self.c_1 = None
        self.c_2 = None

        self.update()

    def update(self):
        assert(self.c.lb >= 0 and self.d.lb >= 0)
        
        if self.c_1 != None:
            self.model.remove(self.c_1)
        if self.c_2 != None:
            self.model.remove(self.c_2)

        p1_a = self.a.lb
        p1_b = self.b.lb
        p1_c = self.c.lb
        p1_d = (p1_a**2 + p1_b**2)/p1_c

        #self.c_1 = self.model.addConstr(0 >= 2*p1_a*(self.a - p1_a) + 2*p1_b*(self.b - p1_b) - p1_d*(self.c - p1_c) - p1_c*(self.d - p1_d), self.prefix+'soc_ub_1'+self.postfix)


        p1_a = self.a.lb
        p1_b = (self.b.ub + self.b.lb)/2
        p1_c = self.c.lb
        p1_d = (p1_a**2 + p1_b**2)/p1_c

        self.c_2 = self.model.addConstr(0 >= 2*p1_a*(self.a - p1_a) + 2*p1_b*(self.b - p1_b) - p1_d*(self.c - p1_c) - p1_c*(self.d - p1_d), self.prefix+'soc_ub_1'+self.postfix)


class LinComplexSOCRelaxScheme(object):
    # relaxation for the equation wr^2 + wi^2 - wfb*wtb <= 0
    def __init__(self, model, wf, wt, wr, wi, vf, vt, td, prefix='', postfix=''):
        self.model = model
        self.prefix = prefix
        self.postfix = postfix

        self.wf = wf
        self.wt = wt
        self.wr = wr
        self.wi = wi

        self.vf = vf
        self.vt = vt
        self.td = td

        self.c_1 = None
        self.c_2 = None

        self.update()

    def update(self):
        assert(self.vf.lb >= 0 and self.vt.lb >= 0)
        
        if self.c_1 != None:
            self.model.remove(self.c_1)
        if self.c_2 != None:
            self.model.remove(self.c_2)

        vf_mp = (self.vf.ub + self.vf.lb)/2
        vt_mp = (self.vt.ub + self.vt.lb)/2
        td_mp = (self.td.ub + self.td.lb)/2

        p_wr = vf_mp*vt_mp*cos(td_mp)
        p_wi = vf_mp*vt_mp*sin(td_mp)
        p_wf = vf_mp**2
        p_wt = vt_mp**2
        
        self.c_1 = self.model.addConstr(0 >= 2*p_wr*(self.wr - p_wr) + 2*p_wi*(self.wi - p_wi) - p_wt*(self.wf - p_wf) - p_wf*(self.wt - p_wt), self.prefix+'soc_ub_1'+self.postfix)
        
        # p_wr = self.vf.ub*self.vt.ub*cos(td_mp)
        # p_wi = self.vf.ub*self.vt.ub*sin(td_mp)
        # p_wf = self.vf.ub**2
        # p_wt = self.vt.ub**2
        
        # self.c_1 = self.model.addConstr(0 >= 2*p_wr*(self.wr - p_wr) + 2*p_wi*(self.wi - p_wi) - p_wt*(self.wf - p_wf) - p_wf*(self.wt - p_wt), self.prefix+'soc_ub_1'+self.postfix)
        
        # p_wr = self.vf.lb*self.vt.lb*cos(td_mp)
        # p_wi = self.vf.lb*self.vt.lb*sin(td_mp)
        # p_wf = self.vf.lb**2
        # p_wt = self.vt.lb**2
        
        # self.c_2 = self.model.addConstr(0 >= 2*p_wr*(self.wr - p_wr) + 2*p_wi*(self.wi - p_wi) - p_wt*(self.wf - p_wf) - p_wf*(self.wt - p_wt), self.prefix+'soc_ub_1'+self.postfix)



class QCModel(object):
    def __init__(self, model, case, output_level=0):
        self.output_level = output_level

        self.case = case
        self.m = model

        self.ref_bus = case.ref_bus().bus_i
        #print('ref bus', ref_bus)

        self.buses = {x.bus_i : x for x in case.bus}

        self.bus_pairs = set()

        self.arcs_from = {}
        self.arcs_to = {}
        self.arcs = {}
        for x in case.branch:
            self.arcs_from[(x.idx, x.fbus, x.tbus)] = x
            self.arcs_to[(x.idx, x.tbus, x.fbus)] = x
            self.bus_pairs.add((x.fbus, x.tbus))

        self.arcs = self.arcs_from.copy()
        self.arcs.update(self.arcs_to)
        #print bus_pairs

        self.bp_line = {}
        self.bp_ad_lb = {}
        self.bp_ad_ub = {}
        for (l,fb,tb), br in self.arcs_from.iteritems():
            bp = (fb, tb)
            
            if bp in self.bp_ad_lb.keys():
                self.bp_ad_lb[bp] = max(br.angmin, self.bp_ad_lb[bp])
            else:
                self.bp_ad_lb[bp] = br.angmin

            if bp in self.bp_ad_ub.keys():
                self.bp_ad_ub[bp] = min(br.angmax, self.bp_ad_ub[bp])
            else:
                self.bp_ad_ub[bp] = br.angmax

            if bp in self.bp_line.keys():
                self.bp_line[bp] = min(l, self.bp_line[bp])
            else:
                self.bp_line[bp] = l


        self.gens = {x.idx : x for x in case.gen}
        self.gencosts = {x.idx : x for x in case.gencost}

        self.bus_gens = {b:[] for b in self.buses.keys()}
        for g,gen in self.gens.iteritems():
            b = gen.bus
            self.bus_gens[b].append(g)

        self.bus_arcs = {b:[] for b in self.buses.keys()}
        for (l,fb,tb) in self.arcs.keys():
            self.bus_arcs[fb].append((l,fb,tb))

        self.gencost_expr = None
        self.const_schemes = []
        self.const_schemes_2 = []
        self.boundchanged = set()
        self.var_schemes = {}
        self.var_schemes_2 = {}
        #print bus_arcs

    def addSchemeVars(self, scheme, svars):
        for var in svars:
            if not var.varName in self.var_schemes.keys():
                self.var_schemes[var.varName] = set()
            self.var_schemes[var.varName].add(scheme)

    def addScheme2Vars(self, scheme, svars):
        for var in svars:
            if not var.varName in self.var_schemes_2.keys():
                self.var_schemes_2[var.varName] = set()
            self.var_schemes_2[var.varName].add(scheme)


    def boundChanged(self, var):
        self.boundchanged.add(var)

    def setObjective(self, expr, sense):
        self.m.setObjective(expr, sense)
        return expr

    def setMinCostObjective(self):
        return self.setObjective(self.gencost_expr, GRB.MINIMIZE)

    def update(self):
        update_time = time.time()

        touched_const_schemes = set()
        touched_const_schemes_2 = set()
        
        if self.output_level >= 5:
            print('vars touched:', len(self.boundchanged))
        for v in self.boundchanged:
            touched_const_schemes |= self.var_schemes[v.varName]
            touched_const_schemes_2 |= self.var_schemes_2[v.varName]
        if self.output_level >= 5:
            print('constr touched:', len(touched_const_schemes), len(touched_const_schemes_2))

        self.m.update()
        for const in touched_const_schemes:
            const.update()
        self.m.update()
        for const in touched_const_schemes_2:
            const.update()
        self.m.update()
        self.boundchanged.clear()

        if self.output_level >= 5:
            print('update time:', time.time() - update_time)

    def update_all(self):
        update_time = time.time()

        self.m.update()
        for const in self.const_schemes:
            const.update()
        self.m.update()
        for const in self.const_schemes_2:
            const.update()
        self.m.update()
        self.boundchanged.clear()

        if self.output_level >= 5:
            print('update time:', time.time() - update_time)

    def build(self):
        m = self.m 
        case = self.case

        time_0 = time.time()
        if self.output_level >= 5:
            print('adding variables')

        self.pg = {i : m.addVar(
            lb=x.Pmin, ub=x.Pmax, 
            vtype=GRB.CONTINUOUS, 
            name='pg_'+str(x.bus)+'_'+str(i)) 
            for i,x in self.gens.iteritems()}
            
        self.qg = {i : m.addVar(
            lb=x.Qmin, ub=x.Qmax, 
            vtype=GRB.CONTINUOUS, 
            name='qg_'+str(x.bus)+'_'+str(i)) 
            for i,x in self.gens.iteritems()}

        p = {}
        q = {}
        for bid, br in self.arcs.iteritems():
            p[bid] = m.addVar(
                lb=-br.rateA, ub=br.rateA, 
                vtype=GRB.CONTINUOUS, 
                name='p_'+str(bid))
            q[bid] = m.addVar(
                lb=-br.rateA, ub=br.rateA, 
                vtype=GRB.CONTINUOUS, 
                name='q_'+str(bid))

        self.v = {i : m.addVar(
            lb=x.Vmin, ub=x.Vmax, 
            vtype=GRB.CONTINUOUS, 
            name='v_'+str(i)) 
            for i,x in self.buses.iteritems()}

        self.t = {i : m.addVar(
            lb=-float('inf'), ub=float('inf'), 
            vtype=GRB.CONTINUOUS, 
            name='t_'+str(i)) 
            for i,x in self.buses.iteritems()}

        w = {i : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='w_'+str(i)) 
            for i,x in self.buses.iteritems()}

        self.td = {(i,j) : m.addVar(
            lb=self.bp_ad_lb[(i,j)], ub=self.bp_ad_ub[(i,j)], 
            vtype=GRB.CONTINUOUS, 
            name='td_'+str((i,j))) 
            for i,j in self.bus_pairs}

        wr = {(i,j) : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='wr_'+str((i,j))) 
            for i,j in self.bus_pairs}

        wi = {(i,j) : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='wi_'+str((i,j))) 
            for i,j in self.bus_pairs}

        vv = {(i,j) : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='vv_'+str((i,j))) 
            for i,j in self.bus_pairs}

        cs = {(i,j) : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='cs_'+str((i,j))) 
            for i,j in self.bus_pairs}

        si = {(i,j) : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='si_'+str((i,j))) 
            for i,j in self.bus_pairs}

        c = {(i,j) : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='c_'+str((i,j))) 
            for i,j in self.bus_pairs}


        if self.output_level >= 5:
            print('model update', time.time() - time_0)
        m.update()

        if self.output_level >= 5:
            print('building fule cost objective', time.time() - time_0)
        self.gencost_expr = 0
        for g,gen in self.gens.iteritems():
            cost = self.gencosts[g]
            self.gencost_expr += cost.c2*self.pg[g]*self.pg[g] + cost.c1*self.pg[g] + cost.c0

        if self.output_level >= 5:
            print('adding model constraints', time.time() - time_0)
        m.addConstr(self.t[self.ref_bus] == 0, 'reference_bus')

        if self.output_level >= 5:
            print('  bus constraints', time.time() - time_0)
        for b,bus in self.buses.iteritems():
            const_s = SquareRelaxScheme(m, self.v[b], w[b])
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (self.v[b], w[b]) )

            p_gen = 0
            q_gen = 0

            for g in self.bus_gens[b]:
                #print b, g, gens[g]
                p_gen += self.pg[g]
                q_gen += self.qg[g]

            p_flow = 0
            q_flow = 0

            for (l,fb,tb) in self.bus_arcs[b]:
                p_flow += p[(l,fb,tb)]
                q_flow += q[(l,fb,tb)]

            #print('bus', bus.Pd, bus.Qd, bus.Gs, bus.Bs)
            m.addConstr(p_flow + bus.Gs*w[b] + bus.Pd == p_gen, 'KCL_P_'+str(b))
            m.addConstr(q_flow - bus.Bs*w[b] + bus.Qd == q_gen, 'KCL_Q_'+str(b))

        del b, g

        if self.output_level >= 5:
            print('  line constraints', time.time() - time_0)
        for (l,fb,tb), br in self.arcs_from.iteritems():
            fbid = (l,fb,tb)
            tbid = (l,tb,fb)
            bp = (fb,tb)
            #print fbid, bp 
            
            gl =  br.r / (br.r**2 + br.x**2)
            bl = -br.x / (br.r**2 + br.x**2)

            tr = br.ratio
            if tr == 0:
                tr = 1.0

            cc = tr*cos(br.angle)
            dd = tr*sin(br.angle)
            tr2 = tr**2

            m.addConstr(p[fbid] == gl/tr2*w[fb] + (-gl*cc + bl*dd)/tr2*wr[bp] + (-bl*cc - gl*dd)/tr2*wi[bp], 'P_'+str(fbid))
            m.addConstr(p[tbid] == gl*w[tb]     + (-gl*cc - bl*dd)/tr2*wr[bp] + (-bl*cc + gl*dd)/tr2*-wi[bp], 'P_'+str(tbid))

            m.addConstr(q[fbid] == -(br.b/2 + bl)/tr2*w[fb] - (-bl*cc - gl*dd)/tr2*wr[bp] + (-gl*cc + bl*dd)/tr2*wi[bp], 'Q_'+str(fbid))
            m.addConstr(q[tbid] == -(br.b/2 + bl)*w[tb]     - (-bl*cc + gl*dd)/tr2*wr[bp] + (-gl*cc - bl*dd)/tr2*-wi[bp], 'Q_'+str(tbid))

            m.addConstr(p[fbid]*p[fbid] + q[fbid]*q[fbid] <= br.rateA**2, 's_'+str(fbid))
            m.addConstr(p[tbid]*p[tbid] + q[tbid]*q[tbid] <= br.rateA**2, 's_'+str(tbid))

        if self.output_level >= 5:
            print('  bus pair constraints', time.time() - time_0)
        for fb,tb in self.bus_pairs:
            bp = (fb,tb)
            postfix = '_'+str(bp)

            m.addConstr(self.td[bp] == self.t[fb] - self.t[tb], 'td_link_'+str(bp))

            const_s = SineRelaxScheme(m, self.td[bp], si[bp], postfix=postfix)
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (self.td[bp], si[bp]) )

            const_s = CosineRelaxScheme(m, self.td[bp], cs[bp], postfix=postfix)
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (self.td[bp], cs[bp]) )

            const_s = ProdRelaxScheme(m, self.v[fb], self.v[tb], vv[bp], prefix='vv_', postfix=postfix)
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (self.v[fb], self.v[tb], vv[bp]) )

            const_s = WPADRelaxScheme(m, wr[bp], wi[bp], self.td[bp], postfix=postfix)
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (wr[bp], wi[bp], self.td[bp]) )

            #useful for testing
            #m.addConstr(wr[bp]*wr[bp] + wi[bp]*wi[bp] <= w[fb]*w[tb], 'w_prod_ub_'+str(bp))

            l = self.bp_line[bp] 
            fbid = (l,fb,tb)
            br = self.arcs_from[fbid]

            gl =  br.r / (br.r**2 + br.x**2)
            bl = -br.x / (br.r**2 + br.x**2)

            tr = br.ratio
            if tr == 0:
                tr = 1.0

            cc = tr*cos(br.angle)
            dd = tr*sin(br.angle)
            tr2 = tr*tr

            const_s = CurrentSqrRelaxScheme(m, self.v[fb], c[bp], br.rateA, tr)
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (self.v[fb], c[bp]) )

            m.addConstr(p[fbid]*p[fbid] + q[fbid]*q[fbid] <= w[fb]/tr2*c[bp], 'l_ub_'+str(bp))
            m.addConstr(c[bp] == (gl*gl + bl*bl)*(w[fb]/tr2 + w[tb] - 2*(cc*wr[bp] + dd*wi[bp])/tr2) - br.b*q[fbid] - ((br.b/2)/tr)**2*w[fb], 'l_link_'+str(bp))

            const_s = LNCScheme(m, w[fb], w[tb], wr[bp], wi[bp], self.v[fb], self.v[tb], self.td[bp], postfix=postfix)
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (w[fb], w[tb], wr[bp], wi[bp], self.v[fb], self.v[tb], self.td[bp]) )

        if self.output_level >= 5:
            print('  model update', time.time() - time_0)
        m.update()

        #needs to be built after update, so the bounds on vv, si, cs, are update to date.
        if self.output_level >= 5:
            print('  bus pair w constraints', time.time() - time_0)
        for fb,tb in self.bus_pairs:
            bp = (fb,tb)
            postfix = '_'+str(bp)

            const_s = ProdRelaxScheme(m, vv[bp], si[bp], wi[bp], prefix='wi_', postfix=postfix)
            self.const_schemes_2.append(const_s)
            self.addScheme2Vars(const_s, (vv[bp], si[bp], wi[bp], self.v[fb], self.v[tb], self.td[bp]) )

            const_s = ProdRelaxScheme(m, vv[bp], cs[bp], wr[bp], prefix='wr_', postfix=postfix)
            self.const_schemes_2.append(const_s)
            self.addScheme2Vars(const_s, (vv[bp], cs[bp], wr[bp], self.v[fb], self.v[tb], self.td[bp]) )

        if self.output_level >= 5:
            print('model update', time.time() - time_0)
        m.update()


    def build_linear(self):
        m = self.m 
        case = self.case

        time_0 = time.time()
        if self.output_level >= 5:
            print('adding variables')

        self.pg = {i : m.addVar(
            lb=x.Pmin, ub=x.Pmax, 
            vtype=GRB.CONTINUOUS, 
            name='pg_'+str(x.bus)+'_'+str(i)) 
            for i,x in self.gens.iteritems()}
            
        self.qg = {i : m.addVar(
            lb=x.Qmin, ub=x.Qmax, 
            vtype=GRB.CONTINUOUS, 
            name='qg_'+str(x.bus)+'_'+str(i)) 
            for i,x in self.gens.iteritems()}

        p = {}
        q = {}
        for bid, br in self.arcs.iteritems():
            p[bid] = m.addVar(
                lb=-br.rateA, ub=br.rateA, 
                vtype=GRB.CONTINUOUS, 
                name='p_'+str(bid))
            q[bid] = m.addVar(
                lb=-br.rateA, ub=br.rateA, 
                vtype=GRB.CONTINUOUS, 
                name='q_'+str(bid))

        self.v = {i : m.addVar(
            lb=x.Vmin, ub=x.Vmax, 
            vtype=GRB.CONTINUOUS, 
            name='v_'+str(i)) 
            for i,x in self.buses.iteritems()}

        self.t = {i : m.addVar(
            lb=-float('inf'), ub=float('inf'), 
            vtype=GRB.CONTINUOUS, 
            name='t_'+str(i)) 
            for i,x in self.buses.iteritems()}

        w = {i : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='w_'+str(i)) 
            for i,x in self.buses.iteritems()}

        self.td = {(i,j) : m.addVar(
            lb=self.bp_ad_lb[(i,j)], ub=self.bp_ad_ub[(i,j)], 
            vtype=GRB.CONTINUOUS, 
            name='td_'+str((i,j))) 
            for i,j in self.bus_pairs}

        wr = {(i,j) : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='wr_'+str((i,j))) 
            for i,j in self.bus_pairs}

        wi = {(i,j) : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='wi_'+str((i,j))) 
            for i,j in self.bus_pairs}

        vv = {(i,j) : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='vv_'+str((i,j))) 
            for i,j in self.bus_pairs}

        cs = {(i,j) : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='cs_'+str((i,j))) 
            for i,j in self.bus_pairs}

        si = {(i,j) : m.addVar(
            vtype=GRB.CONTINUOUS, 
            name='si_'+str((i,j))) 
            for i,j in self.bus_pairs}

        if self.output_level >= 5:
            print('model update', time.time() - time_0)
        m.update()

        if self.output_level >= 5:
            print('building fule cost objective', time.time() - time_0)
        self.gencost_expr = 0
        for g,gen in self.gens.iteritems():
            cost = self.gencosts[g]
            self.gencost_expr += cost.c2*self.pg[g]*self.pg[g] + cost.c1*self.pg[g] + cost.c0

        if self.output_level >= 5:
            print('adding model constraints', time.time() - time_0)
        m.addConstr(self.t[self.ref_bus] == 0, 'reference_bus')

        if self.output_level >= 5:
            print('  bus constraints', time.time() - time_0)
        for b,bus in self.buses.iteritems():
            #const_s = SquareRelaxScheme(m, self.v[b], w[b])
            const_s = LinSquareRelaxScheme(m, self.v[b], w[b])
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (self.v[b], w[b]) )

            p_gen = 0
            q_gen = 0

            for g in self.bus_gens[b]:
                #print b, g, gens[g]
                p_gen += self.pg[g]
                q_gen += self.qg[g]

            p_flow = 0
            q_flow = 0

            for (l,fb,tb) in self.bus_arcs[b]:
                p_flow += p[(l,fb,tb)]
                q_flow += q[(l,fb,tb)]

            #print('bus', bus.Pd, bus.Qd, bus.Gs, bus.Bs)
            m.addConstr(p_flow + bus.Gs*w[b] + bus.Pd == p_gen, 'KCL_P_'+str(b))
            m.addConstr(q_flow - bus.Bs*w[b] + bus.Qd == q_gen, 'KCL_Q_'+str(b))

        del b, g

        if self.output_level >= 5:
            print('  line constraints', time.time() - time_0)
        for (l,fb,tb), br in self.arcs_from.iteritems():
            fbid = (l,fb,tb)
            tbid = (l,tb,fb)
            bp = (fb,tb)
            #print fbid, bp 
            
            gl =  br.r / (br.r**2 + br.x**2)
            bl = -br.x / (br.r**2 + br.x**2)

            tr = br.ratio
            if tr == 0:
                tr = 1.0

            cc = tr*cos(br.angle)
            dd = tr*sin(br.angle)
            tr2 = tr**2

            #print g, b, cc, dd
            m.addConstr(p[fbid] == gl/tr2*w[fb] + (-gl*cc + bl*dd)/tr2*wr[bp] + (-bl*cc - gl*dd)/tr2*wi[bp], 'P_'+str(fbid))
            m.addConstr(p[tbid] == gl*w[tb]     + (-gl*cc - bl*dd)/tr2*wr[bp] + (-bl*cc + gl*dd)/tr2*-wi[bp], 'P_'+str(tbid))

            m.addConstr(q[fbid] == -(br.b/2 + bl)/tr2*w[fb] - (-bl*cc - gl*dd)/tr2*wr[bp] + (-gl*cc + bl*dd)/tr2*wi[bp], 'Q_'+str(fbid))
            m.addConstr(q[tbid] == -(br.b/2 + bl)*w[tb]     - (-bl*cc + gl*dd)/tr2*wr[bp] + (-gl*cc - bl*dd)/tr2*-wi[bp], 'Q_'+str(tbid))

            # m.addConstr(p[fbid]*p[fbid] + q[fbid]*q[fbid] <= br.rateA**2, 's_'+str(fbid))
            # m.addConstr(p[tbid]*p[tbid] + q[tbid]*q[tbid] <= br.rateA**2, 's_'+str(tbid))
            const_s = LinCircleRelaxScheme(m, p[fbid], q[fbid], br.rateA, postfix=str(fbid))
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (p[fbid], q[fbid]) )
            
            const_s = LinCircleRelaxScheme(m, p[tbid], q[tbid], br.rateA, postfix=str(tbid))
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (p[tbid], q[tbid]) )

        if self.output_level >= 5:
            print('  bus pair constraints', time.time() - time_0)
        for fb,tb in self.bus_pairs:
            bp = (fb,tb)
            postfix = '_'+str(bp)

            m.addConstr(self.td[bp] == self.t[fb] - self.t[tb], 'td_link_'+str(bp))

            const_s = SineRelaxScheme(m, self.td[bp], si[bp], postfix=postfix)
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (self.td[bp], si[bp]) )

            #const_s = CosineRelaxScheme(m, self.td[bp], cs[bp], postfix=postfix)
            const_s = LinCosineRelaxScheme(m, self.td[bp], cs[bp], postfix=postfix)
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (self.td[bp], cs[bp]) )

            const_s = ProdRelaxScheme(m, self.v[fb], self.v[tb], vv[bp], prefix='vv_', postfix=postfix)
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (self.v[fb], self.v[tb], vv[bp]) )

            #const_s = WPADRelaxScheme(m, wr[bp], wi[bp], self.td[bp], postfix=postfix)
            #self.const_schemes.append(const_s)
            #self.addSchemeVars(const_s, (wr[bp], wi[bp], self.td[bp]) )


            #useful for testing
            # m.addConstr(wr[bp]*wr[bp] + wi[bp]*wi[bp] <= w[fb]*w[tb], 'w_prod_ub_'+str(bp))

            # const_s = LinComplexSOCRelaxScheme(m, w[fb], w[tb], wr[bp], wi[bp], self.v[fb], self.v[tb], self.td[bp], postfix=postfix)
            # self.const_schemes.append(const_s)
            # self.addSchemeVars(const_s, (w[fb], w[tb], wr[bp], wi[bp], self.v[fb], self.v[tb], self.td[bp]))

            const_s = LNCScheme(m, w[fb], w[tb], wr[bp], wi[bp], self.v[fb], self.v[tb], self.td[bp], postfix=postfix)
            self.const_schemes.append(const_s)
            self.addSchemeVars(const_s, (w[fb], w[tb], wr[bp], wi[bp], self.v[fb], self.v[tb], self.td[bp]) )

        if self.output_level >= 5:
            print('  model update', time.time() - time_0)
        m.update()

        #needs to be built after update, so the bounds on vv, si, cs, are update to date.
        if self.output_level >= 5:
            print('  bus pair w constraints', time.time() - time_0)
        for fb,tb in self.bus_pairs:
            bp = (fb,tb)
            postfix = '_'+str(bp)

            const_s = ProdRelaxScheme(m, vv[bp], si[bp], wi[bp], prefix='wi_', postfix=postfix)
            self.const_schemes_2.append(const_s)
            self.addScheme2Vars(const_s, (vv[bp], si[bp], wi[bp], self.v[fb], self.v[tb], self.td[bp]) )

            const_s = ProdRelaxScheme(m, vv[bp], cs[bp], wr[bp], prefix='wr_', postfix=postfix)
            self.const_schemes_2.append(const_s)
            self.addScheme2Vars(const_s, (vv[bp], cs[bp], wr[bp], self.v[fb], self.v[tb], self.td[bp]) )

            l = self.bp_line[bp] 
            fbid = (l,fb,tb)
            br = self.arcs_from[fbid]

            const_s = CurrentSqrBoundsRelaxScheme(m, w[fb], w[tb], wr[bp], wi[bp], q[fbid], br, postfix=postfix)
            self.const_schemes_2.append(const_s)
            self.addScheme2Vars(const_s, (w[fb], w[tb], wr[bp], wi[bp], q[fbid]) )

            #useful for testing
            #m.addConstr(wr[bp]*wr[bp] + wi[bp]*wi[bp] <= w[fb]*w[tb], 'w_prod_ub_'+str(bp))
            
            #const_s = LinSOCRelaxScheme(m, wr[bp], wi[bp], w[fb], w[tb], postfix=postfix)
            #self.const_schemes_2.append(const_s)
            #self.addScheme2Vars(const_s, (w[fb], w[tb], wr[bp], wi[bp]))

        if self.output_level >= 5:
            print('model update', time.time() - time_0)
        m.update()


