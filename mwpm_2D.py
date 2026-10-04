import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from tqdm import tqdm
import itertools
import math



def create_noise(px, pz, n_qubits):
    list_noise = np.full(n_qubits, 'I')

    for i in range(n_qubits):
        if np.random.random() < px:
            list_noise[i] = 'X'

    return list_noise

#TODO vortex
def apply_noise(list_noise, qc):
    qubit_index = 0

    for noise in list_noise:

        if(noise == 'X'):
            qc.x(qubit_index)
        
        qubit_index += 1
    qc.barrier()


def detection_syndrom(n_ancilla, qc, n_qubits_x, n_qubits_y, simulator, graph_circuit):
    n_data = qc.num_qubits - n_ancilla
    for i in range(n_qubits_x):
        for j in range(n_qubits_y):
            reference = i*(2*n_qubits_x - 1) + j
            ancilla_indice = n_data + i*n_qubits_y + j
            graph_circuit.add_node(
                ancilla_indice,
                type="ancilla"
            )

            qc.cx(reference, ancilla_indice)
            graph_circuit.add_edge(reference, ancilla_indice)
            graph_circuit.nodes[reference]["type"] = "qubit"

            qc.cx(reference + 1, ancilla_indice)
            graph_circuit.add_edge(reference + 1, ancilla_indice)
            graph_circuit.nodes[reference + 1]["type"] = "qubit"

            if j == 0:
                border = ('B', ancilla_indice)

                graph_circuit.add_edge(reference, border)
                graph_circuit.nodes[border]["type"] = "border"


            if j == n_qubits_y - 1:
                border = ('B', ancilla_indice)

                graph_circuit.add_edge(reference + 1, border)
                graph_circuit.nodes[border]["type"] = "border"
                
            if(i > 0):
                index = reference - n_qubits_x + 1
                qc.cx(index, ancilla_indice)
                graph_circuit.add_edge(index, ancilla_indice)
                graph_circuit.nodes[index]["type"] = "qubit"

            if(i < n_qubits_y):
                index = reference + n_qubits_x
                qc.cx(index, ancilla_indice)
                graph_circuit.add_edge(index, ancilla_indice)
                graph_circuit.nodes[index]["type"] = "qubit"

            qc.barrier()


    for i in range(n_ancilla):
        qc.measure(n_data + i, i)


    result = simulator.run(qc, shots=1).result()
    bitstring = list(result.get_counts().keys())[0]
    syndrom = [int(bit) for bit in bitstring[::-1]]

    return syndrom, graph_circuit
    

def make_graph(syndrom, graph, n_qubits_y, n_qubits_x):
    n_x = n_qubits_x - 1 
    nodes = [i for i, v in enumerate(syndrom) if v == 1]
    graph = nx.complete_graph(nodes)
    for i, j in graph.edges:
        ri, ci = divmod(i, n_x)
        rj, cj = divmod(j, n_x)
        graph[i][j]["weight"] = abs(ri - rj) + abs(ci - cj)

    for i in nodes:
        col = i % n_x
        graph.add_edge(i, ("B", i), weight=min(col + 1, n_x - col))

    for u, v in itertools.combinations(nodes, 2):
        graph.add_edge(("B", u), ("B", v), weight=0)

    matching = nx.min_weight_matching(graph)
    return graph, matching


def create_correction(matching, graph_circuit):
    list_correction = np.full(n_qubits, 'I')

    for i, j in matching:
        if(type(i) == tuple and type(j) == tuple):
            continue
        if(type(i) == tuple):
            print(qc.num_qubits)
            i = ('B', i[1] + qc.num_qubits - n_ancilla)
        else :
            i += qc.num_qubits - n_ancilla

        if(type(j) == tuple):
            j = ('B', j[1] + qc.num_qubits - n_ancilla)
        else :
            j += qc.num_qubits - n_ancilla

        path = nx.shortest_path(graph_circuit, i, j)
        print(path)
        path_qubit = [node for node in path if(graph_circuit.nodes[node]["type"] == "qubit")]
        for qubit in path_qubit:
            list_correction[qubit] = 'X'

    return list_correction








n_qubits_x = 5
n_qubits_y = 4
n_ancilla = n_qubits_x * n_qubits_y
qc = QuantumCircuit((n_qubits_x*(n_qubits_y + 1) + n_qubits_y**2) + n_ancilla, n_ancilla)



n_qubits = (n_qubits_x*(n_qubits_y + 1) + n_qubits_y**2)
simulator = AerSimulator()
list_noise = create_noise(0.1 , None, n_qubits)
apply_noise(list_noise, qc)
syndrom, graph_circuit = detection_syndrom(n_ancilla, qc, n_qubits_x, n_qubits_y, simulator, nx.Graph())
print(qc)

graph, matching = make_graph(syndrom , nx.Graph(), n_qubits_y, n_qubits_x)
print(f" matchin : {matching}")
list_correction = create_correction(matching, graph_circuit)
print(list_noise)

print(list_correction)
