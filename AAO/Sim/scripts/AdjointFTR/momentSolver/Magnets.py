class QuadProfile:

    def __init__(self, func, sf=1.0):
        self.func = func
        self.sf = sf

    def GetValue(self, x):
        return (self.func(x) * self.sf)

    def SetScaleFactor(self, sf):
        self.sf = sf

class SolenoidProfile:

    def __init__(self, func, sf=1.0):
        self.func = func
        self.sf = sf

    def GetValue(self, x):
        return ( self.func(x, self.sf) )

    def SetScaleFactor(self, sf):
        self.sf = sf