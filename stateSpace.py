
from copy import deepcopy
from math import sqrt
from base import *

class state:
    def __init__(self, pi:dict) -> None:
        self.m = pi.copy()
    def copy(self):
        s = state(self.m)
        return s
    def actions(self, Q):
        oVs = self.m.values()
        q = len(oVs)
        A = []
        if q == len(Q):
            return A
        for v in range(nHardwareQubit):
            if v not in oVs:
                A.append([q, v])
        return A
    def act(self, a):
        self.m[a[0]] = a[1]
    def isTerminate(self):
        if self.actions() == []:
            return True
        else:
            return False
    def print(self):
        print(self.m)

class tree:
    def __init__(self, s:state) -> None:
        self.s = s
        self.n = 0
        self.v = -1e10
        self.tv = -1e10
        self.children = []
        self.parent = None
    def expand(self, c):
        self.children.append(c)
        c.parent = self

    def print(self, layer=0):
        # if layer > 2:
        #     return
        if self.n == 0:
            return
        for i in range(layer):
            print('   ', end='')
        self.s.print()
        for i in range(layer):
            print('   ', end='')
        print('v:{0} n:{1}'.format(self.v, self.n))
        for c in self.children:
            c.print(layer + 1)

    def depth(self):
        if len(self.children) == 0:
            return 0
        return max([c.depth() for c in self.children]) + 1