import argparse
import numpy as np
import networkx as nx
try:
    import matplotlib.pyplot as plt
except:
    raise


class Kirchhoff():
    def __init__(self):
        # Create new directed graph
        self.graph = nx.DiGraph()

        self.solution = None
        self.SEM = None
        self.eps = 10e-10;

    def load(self, filename):
        # Open file
        f = open(filename, 'r')

        # From first line read vector and value of SEM
        line = f.readline()
        self.SEM = line.split(" ")
        self.SEM[0] = int(self.SEM[0])
        self.SEM[1] = int(self.SEM[1])
        self.SEM[2] = float(self.SEM[2])

        # Edges counter
        edge_num = 0
        node_num = 0

        # Read edges and create graph
        line = f.readline()
        while line and line not in ["", "\n"]:
            splitted = line.split(" ")

            # Edges
            u = int(splitted[0])
            v = int(splitted[1])

            # Resistance beetwen them
            weight = float(splitted[2])

            # Add to graph
            self.graph.add_edge(u, v, weight=weight, edge_num=edge_num)

            edge_num += 1
            line = f.readline()

        for node in self.graph.nodes():
            self.graph.node[node]["node_num"] = node_num
            node_num += 1


        # Close file
        f.close()

    def solve(self):
        # Matrix of coefficients
        a = []

        # Matrix of equation values
        b = []

        # Edges values
        edges_num = nx.get_edge_attributes(self.graph, 'edge_num')

        # Edges amount
        edges_amount = self.graph.number_of_edges()

        # Edges weights
        weights = nx.get_edge_attributes(self.graph, 'weight')

        # First kirchoff law
        for node in self.graph:
            a.append([0]*edges_amount)
            b.append([0])

            # In
            for non_neigh in self.graph.predecessors(node):
                a[-1][edges_num[(non_neigh, node)]] = 1

            # Out
            for neigh in self.graph.successors(node):
                edge = (node, neigh)
                a[-1][edges_num[edge]] = -1

        # Second kirchoff law
        for cycle in nx.cycle_basis(self.graph.to_undirected()):
            # Add last edge like first to close circuit
            cycle.append(cycle[0])

            a.append([0]*edges_amount)

            # If sem is in cycle Sum(I*R) = U otherwise Sum(I*R) = 0
            if (self.SEM[0] in cycle) and (self.SEM[1] in cycle):
                b.append([self.SEM[2]])
            else:
                b.append([0])

            # Go over cycle edges
            for i in range(0, len(cycle)-1):
                # Get edges order which is in graph
                edge_in_graph = edge_in_cycle = (cycle[i], cycle[i+1])

                # Check if direction is like in graph
                if edge_in_cycle not in edges_num:
                    edge_in_graph = (cycle[i+1], cycle[i])

                # Write R to coefficients matrix
                if edge_in_graph == edge_in_cycle:
                    a[-1][edges_num[edge_in_graph]] = (weights[edge_in_graph])
                else:
                    a[-1][edges_num[edge_in_graph]] = -(weights[edge_in_graph])

        # To numpy array
        a = np.array(a)
        b = np.array(b)

        # Solution
        self.solution = (np.linalg.lstsq(a, b))[0]

    def save_i_in_graph(self):
        # Get edges
        edges = self.graph.edges(data=True)

        # Save I in graph
        for edge in edges:
            edge[2]["I"] = self.solution[edge[2]["edge_num"]][0]

    def set_directed_to_flow(self):
        # Get edges
        edges = list(self.graph.edges(data=True))

        # Change flow vector if value is < 0
        for edge in edges:
            # Get edge attrs
            attrs = edge[2]

            if self.solution[attrs["edge_num"]] < 0:
                # Change flow vector
                self.solution[attrs["edge_num"]] *= -1

                # Remove wrong
                self.graph.remove_edge(edge[0], edge[1])

                # Add right with attributes
                self.graph.add_edge(edge[1], edge[0], weight=attrs["weight"], edges_num=attrs["edge_num"], I=-attrs["I"])

    def valid(self):
        non_valid_counter = 0

        # Edges values
        edges_currents = nx.get_edge_attributes(self.graph, 'I')

        for node in self.graph.nodes():
            sum = 0

            for non_neigh in self.graph.predecessors(node):
                sum += edges_currents[(non_neigh, node)]

            # Out
            for neigh in self.graph.successors(node):
                sum -= edges_currents[(node, neigh)]

            if sum > self.eps:
                non_valid_counter += 1

        self.valid = (non_valid_counter == 0)

    def draw(self, filename, spectral):
        # Extract I values
        edge_labels = dict([((u,v,), round(d['I'], 2)) for u, v, d in self.graph.edges(data=True)])

        # Add layout to grapg
        pos=nx.spectral_layout(self.graph) if spectral else nx.circular_layout(self.graph)

        # Create heatmap values
        edges_colors = [abs(round(d['I'])) for u, v, d in self.graph.edges(data=True)]

        # Colors nodes with special color for SEM nodes
        colors = []
        for node in self.graph:
            if node in self.SEM[0:2]:
                colors.append("#ff1a1a")
            else:
                colors.append("#A0CBE2")

        # Set labels on graph
        nx.draw_networkx_edge_labels(self.graph, pos, edge_labels=edge_labels)

        nx.draw(self.graph, pos, node_color=colors, width=1, node_size=750, edge_color=edges_colors, with_labels=True, edge_cmap=plt.cm.Reds)
        plt.savefig(filename)

        # Display info if calculations are valid
        plt.figtext(0.01, 0.01, "Valid: {}\nEpsilon: {}".format(self.valid, self.eps), fontsize=10, color="red")

        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="In data file", type=str, required=True);
    parser.add_argument("--output", help="Output image file", type=str, required=True);
    parser.add_argument("-d", help="Set direction to flow", action='store_true');
    parser.add_argument("-s", help="Set spectral layuot", action='store_true');
    args = parser.parse_args()

    kirchoff = Kirchhoff()
    kirchoff.load(args.input)
    kirchoff.solve()
    kirchoff.save_i_in_graph()

    if args.d:
        kirchoff.set_directed_to_flow()

    kirchoff.valid()

    kirchoff.draw(args.output, args.s)
