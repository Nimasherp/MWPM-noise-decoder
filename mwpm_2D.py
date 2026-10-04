import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from tqdm import tqdm
import itertools
import math


## 1. Returns a list of random noise, for now only bit-flips
def create_noise(px, pz, n_qubits):
    list_noise = np.full(n_qubits, 'I')

    for i in range(n_qubits):
        if np.random.random() < px:
            list_noise[i] = 'X'

    return list_noise

#TODO vortex
## 2. Applies the noise list generated from create_noise 
def apply_noise(list_noise, qc):
    qubit_index = 0

    for noise in list_noise:

        if(noise == 'X'):
            qc.x(qubit_index)
        
        qubit_index += 1
    qc.barrier()

## 3. Detects syndromes and builds the graph-lattice
def detection_syndrom(n_ancilla, qc, n_qubits_x, n_qubits_y, simulator, graph_circuit):
    n_data = qc.num_qubits - n_ancilla

    for i in range(n_qubits_y + 1):
        for j in range(n_qubits_x - 1):

            reference = i*(2*n_qubits_x - 1) + j # indice of the reference qubit
            ancilla_indice = n_data + i*(n_qubits_x - 1) + j
            graph_circuit.add_node(ancilla_indice, type="ancilla")

            qc.cx(reference, ancilla_indice)
            graph_circuit.add_edge(reference, ancilla_indice)
            graph_circuit.nodes[reference]["type"] = "qubit"

            qc.cx(reference + 1, ancilla_indice)
            graph_circuit.add_edge(reference + 1, ancilla_indice)
            graph_circuit.nodes[reference + 1]["type"] = "qubit"

            # Boundary conditions
            if j == 0:
                border = ('B', ancilla_indice)

                graph_circuit.add_edge(reference, border)
                graph_circuit.nodes[border]["type"] = "border"


            if j == n_qubits_x - 2:
                border = ('B', ancilla_indice)

                graph_circuit.add_edge(reference + 1, border)
                graph_circuit.nodes[border]["type"] = "border"

            # Boundary conditions
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


    # Only measures the ancilla
    for i in range(n_ancilla):
        qc.measure(n_data + i, i)


    result = simulator.run(qc, shots=1).result()
    bitstring = list(result.get_counts().keys())[0]
    syndrom = [int(bit) for bit in bitstring[::-1]] # Qiskit uses an upsidedown measurements

    return syndrom, graph_circuit
    
## 4. Decodes error tracking via standard distance weights
def mwpm(syndrom, n_qubits_y, n_qubits_x):
    graph = nx.Graph
    n_x = n_qubits_x - 1 

    nodes = [i for i, v in enumerate(syndrom) if v == 1]
    graph = nx.complete_graph(nodes)

    for i, j in graph.edges: 
        # Values the weight for each syndrom based on the number of qubits separating them
        ri, ci = divmod(i, n_x)
        rj, cj = divmod(j, n_x)
        graph[i][j]["weight"] = abs(ri - rj) + abs(ci - cj)

    for i in nodes:
        # Connects every syndrom to it's border
        col = i % n_x
        graph.add_edge(i, ("B", i), weight=min(col + 1, n_x - col))

    for u, v in itertools.combinations(nodes, 2):
        # Connects every virtual border to eachother with weight = 0
        graph.add_edge(("B", u), ("B", v), weight=0)

    matching = nx.min_weight_matching(graph)
    return matching

## 5. Resets the ancilla to make new measurements after corrections
def reset_ancilla(n_ancilla, qc, n_qubits):
    for i in range(n_ancilla):
        qc.reset(n_qubits + i)
    return qc

## 6. Create a correction list, to apply to the supposed noised qubits
def create_correction(matching, graph_circuit, n_qubits, qc, n_ancilla):
    list_correction = np.full(n_qubits, 'I')
    for i, j in matching:
        if(type(i) == tuple and type(j) == tuple):
            continue
        if(type(i) == tuple):
            i = ('B', i[1] + qc.num_qubits - n_ancilla)
        else :
            i += qc.num_qubits - n_ancilla

        if(type(j) == tuple):
            j = ('B', j[1] + qc.num_qubits - n_ancilla)
        else :
            j += qc.num_qubits - n_ancilla

        path = nx.shortest_path(graph_circuit, i, j)
        path_qubit = [node for node in path if(graph_circuit.nodes[node]["type"] == "qubit")]
        for qubit in path_qubit:
            list_correction[qubit] = 'X'
    return list_correction




def one_experiment(n_qubits_x, n_qubits_y, px, simulator):
    n_ancilla = (n_qubits_x - 1) * (n_qubits_y + 1)
    n_qubits = n_qubits_x*(n_qubits_y + 1) + n_qubits_y*(n_qubits_x - 1)

    qc = QuantumCircuit(n_qubits + n_ancilla, n_ancilla)
    
    list_noise = create_noise(px, None, n_qubits)
    apply_noise(list_noise, qc)

    syndrom, graph_circuit = detection_syndrom(n_ancilla, qc, n_qubits_x, n_qubits_y, simulator, nx.Graph())
    matching = mwpm(syndrom, n_qubits_y, n_qubits_x)
    print(f"syndroms before correction : {syndrom}")
    # visualize_graph(graph_circuit, n_qubits, n_ancilla, n_qubits_x, n_qubits_y)

    list_correction = create_correction(matching, graph_circuit, n_qubits, qc, n_ancilla)
    reset_ancilla(n_ancilla, qc, n_qubits)
    
    apply_noise(list_correction, qc)
    syndrom_corrected, G = detection_syndrom(n_ancilla, qc, n_qubits_x, n_qubits_y, simulator, nx.Graph())
    print(f"syndroms after correction : {syndrom_corrected}")

 


one_experiment(6, 2, 0.1, AerSimulator())

