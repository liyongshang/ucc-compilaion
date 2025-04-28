nHardwareQubit = 53
mChallengeTimes = nHardwareQubit * 9
from graph import load_graph, pGraph
Graph = load_graph('sycamore')
from sys import maxsize
import steiner_forest
steiner_forest.set_graph_data(Graph.adj, Graph.C, nHardwareQubit)

# nHardwareQubit = 100
# mChallengeTimes = nHardwareQubit
# from graph import load_graph, pGraph
# Graph = load_graph('grid1010')
# from sys import maxsize