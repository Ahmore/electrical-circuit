import argparse

import Kirchoff
import numpy as np
import networkx as nx
try:
    import matplotlib.pyplot as plt
except:
    raise


class NodePotential(Kirchoff.Kirchhoff):
    def solve(self):
        # Matrix of coefficients
        a = []

        # Matrix of equation values
        b = []

        # Edges values
        edges_num = nx.get_edge_attributes(self.graph, 'edge_num')

        # Edges weights
        weights = nx.get_edge_attributes(self.graph, 'weight')

        nodes_amount = self.graph.__len__()

        # First kirchoff law
        for node in self.graph:
            a.append([0]*nodes_amount)
            b.append([0])

            # SEM_u as grounded
            if node == self.SEM[0]:
                a[-1][self.graph.node[node]["node_num"]] = 1

            # SEM_u set to U
            elif node == self.SEM[1]:
                a[-1][self.graph.node[node]["node_num"]] = 1
                b[-1] = [self.SEM[2]]

            else:
                for neigh in self.graph.to_undirected().neighbors(node):
                    edge = (node, neigh)
                    if edge not in edges_num:
                        edge = (neigh, node)

                    a[-1][self.graph.node[node]["node_num"]] += 1/weights[edge]
                    a[-1][self.graph.node[neigh]["node_num"]] -= 1/weights[edge]


        # To numpy array
        a = np.array(a)
        b = np.array(b)

        # Solution
        potentials = (np.linalg.lstsq(a, b))[0]

        currents = []

        for i in range(0, self.graph.number_of_edges()):
            currents.append([0])

        for edge in self.graph.edges():
            # Count I on SEM branch from OHMs law
            if edge != (self.SEM[0], self.SEM[1]):
                (u, v) = edge
                vu = potentials[self.graph.node[u]["node_num"]][0]
                vv = potentials[self.graph.node[v]["node_num"]][0]
                currents[edges_num[edge]] = [(vv-vu)/weights[edge]]

        # Count I on SEM branch from 1st Kirchoffs law
        edge = (self.SEM[0], self.SEM[1])

        # In
        for non_neigh in self.graph.predecessors(self.SEM[0]):
            currents[edges_num[edge]][0] += currents[edges_num[(non_neigh, self.SEM[0])]][0]

        # Out
        for neigh in self.graph.successors(self.SEM[0]):
            if edge != (self.SEM[0], neigh):
                currents[edges_num[edge]][0] -= currents[edges_num[(self.SEM[0], neigh)]][0]

        self.solution = np.array(currents)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="In data file", type=str, required=True);
    parser.add_argument("--output", help="Output image file", type=str, required=True);
    parser.add_argument("-d", help="Set direction to flow", action='store_true');
    parser.add_argument("-s", help="Set spectral layuot", action='store_true');
    args = parser.parse_args()

    node_potential = NodePotential()
    node_potential.load(args.input)
    node_potential.solve()
    node_potential.save_i_in_graph()

    if args.d:
        node_potential.set_directed_to_flow()

    node_potential.valid()

    node_potential.draw(args.output, args.s)