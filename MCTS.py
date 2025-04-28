from stateSpace import *
from evaluation import evalPaulis
from tools import rearrange
from sys import maxsize
from time import time
import steiner_forest
import math
class MCTS:
    def __init__(self, Q, eval, _simulation, femi_ops, center):
        reQ, order = rearrange(Q)
        steiner_forest.set_qubit_data(Q, order)
        self.order = order
        self.initState = state({})
        self.initState.act([order[0], center])
        self.tree = tree(self.initState)
        self.eval = eval
        self._simulation = _simulation
        self.Q = Q
        self.B = []
        self.femi_ops = femi_ops
        self.debug = False
        for i in range(len(self.Q[0])):
            b = []
            for j in range(len(self.Q)):
                if self.Q[j][i] == 1:
                    b.append(j)
            self.B.append(b.copy())

    def setChallengeTimes(self, ct):
        self.mChallengeTimes = ct

    def selection(self, t:tree):
        def policy(c:tree):
            alpha = 10
            if c.n == 0:
                return 1.0e20
            return c.tv / c.n + alpha * sqrt(math.log(c.parent.n + 1) / c.n)# + c.v
        t.n += 1
        if len(t.children) == 0:
            return t
        mp = -1e20
        for c in t.children:
            if c != self.champion and policy(c) > mp:
                mc = c
                mp = policy(c)
        return self.selection(mc)
    
    def expansion(self, leaf:tree):
        s = leaf.s
        nmq = len(s.m)
        if nmq >= len(self.Q):
            return
        q = self.order[nmq]
        candV = []
        visited = list(s.m.values())
        for v in visited:
            for _v in Graph.adj[v]:
                if _v not in visited:
                    # visited.append(_v)
                    candV.append(_v)
        # for v in range(nHardwareQubit):
        #     if v not in s.m.values():
        for v in candV:
                newS = s.copy()
                newS.act([q, v])
                leaf.expand(tree(newS))

    def simulation(self, leaf:tree):
        t0 = time()
        m = self._simulation(leaf.s, self.Q, self.order)
        # print('s1: {}'.format(time() - t0))
        t0 = time()
        leaf.v = self.eval(self.femi_ops, m)
        leaf.tv = leaf.v
        # print('s2: {}'.format(time() - t0))
        # input()
        # if leaf.s.m == {0: 10, 1: 6, 2: 9, 3: 11, 4: 13, 5: 14, 6: 12}:
        #     print('test mapping: ', m)
        #     # print(m)
        #     showMappingOnGrid(m, 4, 4)
        #     # exit()
        # if leaf.s.m == {0: 10, 1: 6, 2: 9, 3: 11, 4: 13, 5: 12, 6: 14}:
        #     print('test mapping: ', m)
        #     # print(m)
        #     showMappingOnGrid(m, 4, 4)
        #     # exit()

    def heuristicMapping(self):
        m = self._simulation(self.initState, self.Q, self.order)
        return m

    def backpropagation(self, leaf, simv):
        if leaf == None:
            return
        if len(leaf.children) > 0:
            leaf.v = max(leaf.v, max([c.v for c in leaf.children]))
        leaf.tv += simv
        if leaf.parent != None:
            self.backpropagation(leaf.parent, simv)

    def decision(self):
        if self.debug:
            for c in self.tree.children:
                print(c.v, ':', c.n) # , ' ', end=''
            print('depth:{0}'.format(self.tree.depth()))
            print()
        self.tree = self.champion
        self.tree.parent = None
        # print('decision')
        # self.tree.s.print()
        if self.debug:
            print('v:{0} n:{1} av:{2}'.format(self.tree.v, self.tree.n, self.tree.tv / (self.tree.n + 1.0e-20)))
            input()

    def start(self):
        numLoop = 0
        numUpdateChamp = 0
        simulationTime = 0
        self.simulation(self.tree)
        while(len(self.tree.s.m) < len(self.Q)):
            # print(self.tree.v)
            if len(self.tree.children) == 0:
                self.expansion(self.tree)
                self.champion = None
            if len(self.tree.children) < 2:
                self.champion = self.tree.children[0]
                self.decision()
                continue
            mv = -1e20
            for c in self.tree.children:
                if c.v > mv:
                    mv = c.v
                    self.champion = c
            while(True):
                ct = 0
                while(True):
                    numLoop += 1
                    leaf = self.selection(self.tree)
                    visitedNum = [c.n for c in self.tree.children]
                    totalVisited = sum(visitedNum)
                    self.expansion(leaf)
                    t0 = time()
                    self.simulation(leaf)
                    simulationTime += time() - t0
                    # print('end simulation')
                    self.backpropagation(leaf.parent, leaf.v)
                    ct += 1
                    if leaf.v > self.champion.v:
                        numUpdateChamp += 1
                        break
                    elif ct >= self.mChallengeTimes:
                        break
                    # self.tree.print()
                    # input()
                if ct >= self.mChallengeTimes:
                    break
                newChampion = None
                mv = -1e20
                for c in self.tree.children:
                    if c.v > mv:
                        mv = c.v
                        newChampion = c
                self.champion = newChampion

            visitedNum = [c.n for c in self.tree.children]
            totalVisited = sum(visitedNum)
            if self.debug:
                print('total visited: ', totalVisited)

            self.decision()
        # print(self.tree.v)
        # print(self.eval(self.femi_ops, self.tree.s.m))
        # print('num loop per step: {}\nnum champions per step: {}\nsimulation time: {}'\
            #   .format(numLoop/len(self.Q), numUpdateChamp/len(self.Q), simulationTime))
        self.compileCost = [numLoop, numUpdateChamp]
        