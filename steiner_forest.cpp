#include <vector>
#include <queue>
#include <unordered_set>
#include <tuple>
#include <limits>
#include <unordered_map>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>

namespace py = pybind11;

using namespace std;

// Global graph data
vector<vector<int>> global_G;
vector<vector<int>> global_C;
vector<vector<int>> global_Q;
vector<int> global_order;
int global_nHardwareQubit;
vector<vector<int>> interat(100, vector<int>(100,0));

void set_graph_data(const vector<vector<int>>& G, const vector<vector<int>>& C, int nHardwareQubit) {
    global_G = G;
    global_C = C;
    global_nHardwareQubit = nHardwareQubit;
}

void set_qubit_data(const vector<vector<int>>& Q, const vector<int>& order) {
    global_Q = Q;
    global_order = order;
    vector<double> block_size(Q[0].size(), 0);
    for (int i=0; i<Q.size(); i++){
        for (int j=0; j<Q[0].size(); j++){
            block_size[j] += Q[i][j]
        }
    }
    for (int i=0; i<Q.size(); i++){
        for (int j=0; j<Q.size(); j++){
            for (int k=0; k<Q[i].size(); k++){
                interat[i][j] += Q[i][k] * Q[j][k] / (block_size[k] * block_size[k]);
            }
        }
    }
}

vector<pair<int, int>> SteinerBalanceForest(const vector<int>& roots, const vector<int>& targets) {
    const vector<vector<int>>& G = global_G;
    const vector<vector<int>>& C = global_C;

    vector<pair<int, int>> edges;
    unordered_set<int> addedV(roots.begin(), roots.end());
    vector<int> nodeDegree(G.size(), 0);
    vector<int> nodeDepth(G.size(), 0);
    vector<int> remainTargets = targets;
    vector<pair<int, int>> Ec;

    for (int v : roots) nodeDepth[v] = 1;

    for (int v : roots) {
        for (int _v : G[v]) {
            if (addedV.find(_v) == addedV.end()) {
                Ec.emplace_back(v, _v);
            }
        }
    }

    while (!remainTargets.empty()) {
        pair<int, int> _e;
        int md = 1e9, mp = 1e9;
        int originTreeDepth = *max_element(nodeDepth.begin(), nodeDepth.end());
        int originDegree = *max_element(nodeDegree.begin(), nodeDegree.end());

        for (const auto& e : Ec) {
            int dcost = 0;
            vector<int> tempV(addedV.begin(), addedV.end());
            tempV.push_back(e.second);

            for (int v : targets) {
                if (addedV.find(v) == addedV.end()) {
                    int minDist = 1e9;
                    for (int _v : tempV) {
                        minDist = min(minDist, C[v][_v]);
                    }
                    dcost += minDist;
                }
            }

            int pcost = max({originTreeDepth, nodeDepth[e.first] + 1, originDegree, nodeDegree[e.first] + 1});

            if (dcost < md || (dcost == md && pcost < mp)) {
                _e = e;
                md = dcost;
                mp = pcost;
            }
        }

        edges.push_back(_e);
        addedV.insert(_e.second);
        if (find(remainTargets.begin(), remainTargets.end(), _e.second) != remainTargets.end()) {
            remainTargets.erase(remove(remainTargets.begin(), remainTargets.end(), _e.second), remainTargets.end());
        }
        nodeDegree[_e.first]++;
        nodeDepth[_e.second] = nodeDepth[_e.first] + 1;
        Ec.erase(remove(Ec.begin(), Ec.end(), _e), Ec.end());

        for (int adj : G[_e.second]) {
            if (addedV.find(adj) == addedV.end()) {
                Ec.emplace_back(_e.second, adj);
            }
        }
    }

    return edges;
}

int secondCost(const vector<int>& pobs) {
    const vector<vector<int>>& G = global_G;
    const vector<vector<int>>& C = global_C;

    if (pobs.size() == 2) {
        return C[pobs[0]][pobs[1]] - 1;
    }
    int C2 = numeric_limits<int>::max();

    for (size_t i = 0; i < G.size(); ++i) {
        if (G[i].size() < 3) continue;

        vector<int> tempPobs;
        for (int j : pobs) {
            if (find(G[i].begin(), G[i].end(), j) == G[i].end() && j != i) {
                tempPobs.push_back(j);
            }
        }

        int c = 0;
        int v0;
        if (tempPobs.empty() || find(pobs.begin(), pobs.end(), i) != pobs.end()) {
            v0 = *min_element(pobs.begin(), pobs.end(), [&C, i](int x, int y) { return C[i][x] < C[i][y]; });
        } else {
            v0 = *min_element(tempPobs.begin(), tempPobs.end(), [&C, i](int x, int y) { return C[i][x] < C[i][y]; });
        }
        c += C[i][v0];

        vector<int> candPos = G[i];
        vector<int> remainPobs;

        for (int v : pobs) {
            if (v == v0) continue;
            if (find(candPos.begin(), candPos.end(), v) != candPos.end()) {
                candPos.erase(remove(candPos.begin(), candPos.end(), v), candPos.end());
            } else {
                remainPobs.push_back(v);
            }
        }

        for (int v : remainPobs) {
            auto pos = min_element(candPos.begin(), candPos.end(), [&C, v](int x, int y) { return C[v][x] < C[v][y]; });
            if (pos != candPos.end()) {
                c += C[v][*pos];
                candPos.erase(pos);
            }
        }

        C2 = min(C2, c);
    }
    return C2;
}

unordered_map<int, int> simulation(
    unordered_map<int, int> pi// ,
    // const vector<vector<int>>& Q,
    // const vector<int>& order,
    // int nHardwareQubit
) {
    const vector<vector<int>>& C = global_C;
    const vector<vector<int>>& G = global_G;
    const int& nHardwareQubit = global_nHardwareQubit;
    const vector<vector<int>>& Q = global_Q;
    const vector<int>& order = global_order;
    unordered_set<int> candPos;

    unordered_map<int, int> m = pi;
    size_t nMappedQubits = m.size();

    for (const auto& kv: m){
        int v = kv.second;
        for (const auto& _v: G[v]){
            candPos.insert(_v);
        }
    }

    for (size_t idx = nMappedQubits; idx < order.size(); ++idx) {
        int i = order[idx];
        // const vector<int>& q = Q[i];
        int md = numeric_limits<int>::max();
        int bestVertex = -1;

        for (int vertex = 0; vertex < nHardwareQubit; ++vertex) {
            if (find_if(m.begin(), m.end(), [vertex](const pair<int, int>& p) { return p.second == vertex; }) != m.end()) {
                continue;
            }
            int d = 0;
            for (const auto& kv : m) {
                int k = kv.first;
                int mapped = kv.second;
                // const vector<int>& _q = Q[k];
                // for (size_t j = 0; j < q.size(); ++j) {
                //     d += q[j] * _q[j] * C[vertex][mapped];
                // }
                d += -1 * interat[i][k] / C[vertex][mapped];
            }
            if (d < md) {
                md = d;
                bestVertex = vertex;
            }
        }
        m[i] = bestVertex;
        // candPos.erase(bestVertex);
        // for (const auto& _v: G[bestVertex]){
        //     candPos.insert(_v);
        // }
    }
    return m;
}

PYBIND11_MODULE(steiner_forest, m) {
    m.def("set_graph_data", &set_graph_data, "Initialize G and C graph data");
    m.def("set_qubit_data", &set_qubit_data, "Initialize Q and order qubit data");
    m.def("SteinerBalanceForest", &SteinerBalanceForest, "Constructs a Steiner Balance Forest");
    m.def("secondCost", &secondCost, "Computes the second cost function");
    m.def("simulation", &simulation, "Performs hardware mapping simulation");
}
