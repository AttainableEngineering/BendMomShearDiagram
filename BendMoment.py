import matplotlib.pyplot as pl


class SingularityEqn:
    def __init__(self):
        self.L = 0.0
        self.resolution = 1


class Moment(SingularityEqn):
    def __init__(self, mag, a):
        SingularityEqn.__init__(self)
        # Magnitude and turn on point of vector
        self.m = mag
        self.a = a
        self.x = 0
        [self.M, self.V] = self.iterate()

    def construct(self, m, a):
        meqn = m*(self.x-a)**0
        veqn = 0
        eqns = [meqn, veqn]
        return eqns
    

    def iterate(self):

        # Get and slice length into acceptable resolution
        L = int(SingularityEqn.L)
        res = self.resolution
        r = L*res
        jj = [ii for ii in range(0, r+1)]
        M = []
        V = []
        for points in jj:
            self.x = points/res
            eqns = self.construct(self.m, self.a)
            M.append(eqns[0])
            V.append(eqns[1])
        out = [M, V]
        return out


class PointLoad(SingularityEqn):
    def __init__(self, force, a):
        SingularityEqn.__init__(self)

        self.f = force
        self.a = a
        self.x = 0 
        [self.M, self.V] = self.iterate()

    def construct(self, f, a):
        if self.x >= a:
            meqn = f*(self.x-a)**1
            veqn = f*(self.x-a)**0
            eqns = [meqn, veqn]
        else:
            eqns = [0, 0]
        return eqns
    
    def iterate(self):

        # Get and slice length into acceptable resolution
        L = int(SingularityEqn.L)
        res = self.resolution
        r = L*res
        jj = [ii for ii in range(0, r+1)]
        M = []
        V = []
        for points in jj:
            self.x = points/res
            eqns = self.construct(self.f, self.a)
            M.append(eqns[0])
            V.append(eqns[1])
        out = [M, V]
        return out


class DistributedLoad(SingularityEqn):
    def __init__(self, distload, a, b):
        SingularityEqn.__init__(self)

        self.p = distload
        self.a = a
        self.b = b
        self.x = 0
        [self.M, self.V] = self.iterate()

    def construct(self, p, a, b):
        if self.x >= a and self.x <= b:
            meqn = p/2*(self.x - a)**2
            veqn = p*(self.x - a)
            eqns = [meqn, veqn]
        elif self.x > a and self.x > b:
            meqn = p/2*(self.x - a)**2 - p/2*(self.x - b)**2
            veqn = p*(self.x - a) - p*(self.x - b)
            eqns = [meqn, veqn]
        else:
            eqns = [0,0]
        return eqns
    
    def iterate(self):
        # Get and slice length into acceptable resolution
        L = int(SingularityEqn.L)
        res = self.resolution
        r = L*res
        jj = [ii for ii in range(0, r+1)]
        M = []
        V = []
        for points in jj:
            self.x = points/res
            eqns = self.construct(self.p, self.a, self.b)
            M.append(eqns[0])
            V.append(eqns[1])
        out = [M, V]
        return out

def AssembleEqns(eqns):
    # Combine moments and shear force contributions
    m = []
    v = []
    for eqn in eqns:
        m.append(eqn[0])
        v.append(eqn[1])
    
    # Rearrange and add contributions
    mzip = zip(*m)
    vzip = zip(*v)
    M = [sum(ii) for ii in mzip]
    V = [sum(jj) for jj in vzip]

    eqns = [M, V]
    return eqns

def Plot(eqns, span):
    M = eqns[0]
    V = eqns[1]
    X = span
    pl.figure()
    pl.plot(X, M)
    pl.figure()
    pl.plot(X, V)
    pl.show()


def Main():
    lg = input("Length of bar:")
    SingularityEqn.L = lg

    # a = Moment(-1.0, 10)
    # b = PointLoad(10.0, 50)
    c = DistributedLoad(4, 2, 6)

    # eq = [[a.M, a.V],[b.M, b.V],[c.M, c.V]]
    eq = [[c.M, c.V]]


    r = int(SingularityEqn().resolution)
    inc = int(SingularityEqn.L*r)

    span = [ii for ii in range(0, inc+1)]
    for num in span:
        num = num/r

    eqns = AssembleEqns(eq)

    Plot(eqns, span)






if __name__ == "__main__":
    Main()