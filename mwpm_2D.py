import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from tqdm import tqdm
import math



def create_noise(px, pz, n_qubits_x, n_qubits_y):
    list_noise = np.full((n_qubits_x, n_qubits_y), 'I')

    for i in range(n_qubits_x):
        for j in range(n_qubits_y):
            if np.random.random() < px:
                list_noise[i, j] = 'X'

    return list_noise

def apply_noise(list_noise, qc):
    qubit_index = 0

    for list_noise_x in list_noise:

        for noise in list_noise_x:

            if(noise == 'X'):
                qc.x(qubit_index)
            
            qubit_index += 1

# def detection_syndrom(n_ancilla, qc, n_qubits_x, n_qubits_y, simulator):
#     for i in range(n_qubits_x):
#         for j in range(n_qubits_y):
#             if(i != 0):
#                 reference = i*n_qubits_x*2 - 1 + j
#             else : 
#                 reference = i + j
#             ancilla_indice = n_qubits_x*n_qubits_y + i*n_qubits_x + j
#             print(reference)
#             qc.cx(reference, ancilla_indice)
#             qc.cx(reference + 1, ancilla_indice)
#             if(i > 0):
#                 qc.cx(reference - n_qubits_x - 1, ancilla_indice)
#             if(i < n_qubits_y - 1):
#                 qc.cx(reference + n_qubits_x, ancilla_indice)  

#         # qc.cx()
#         # if((i + 1)%n_qubits_x == 0):
#         #     ligne += 1

#         # qc.cx(i + ligne, n_qubits_y*n_qubits_x + i)
#         # qc.cx(i + ligne + 1, n_qubits_y*n_qubits_x + i)
#         # qc.cx(i + ligne + n_qubits_x, n_qubits_y*n_qubits_x + i)
#         # qc.cx(i + ligne + n_qubits_x + 1, n_qubits_y*n_qubits_x + i)
#     for i in range(n_ancilla):
#         qc.measure(n_qubits_x*n_qubits_y + i, i)

def detection_syndrom(n_ancilla, qc, n_qubits_x, n_qubits_y, simulator):
    n_data = qc.num_qubits - n_ancilla
    for i in range(n_qubits_x):
        for j in range(n_qubits_y):
            reference = i*(2*n_qubits_y + 1) + j
            ancilla_indice = n_data + i*n_qubits_y + j
            qc.cx(reference, ancilla_indice)
            qc.cx(reference + 1, ancilla_indice)
            if(i > 0):
                qc.cx(reference - n_qubits_y, ancilla_indice)
            if(i < n_qubits_x - 1):
                qc.cx(reference + n_qubits_y + 1, ancilla_indice)
    for i in range(n_ancilla):
        qc.measure(n_data + i, i)

    result = simulator.run(qc, shots=1).result()
    bitstring = list(result.get_counts().keys())[0]
    syndrom = [int(bit) for bit in bitstring[::-1]]

    return syndrom
    result = simulator.run(qc, shots=1).result()
    bitstring = list(result.get_counts().keys())[0]
    syndrom = [int(bit) for bit in bitstring[::-1]]

    return syndrom

def make_graph(syndrom, graph, n_qubits_y, n_qubits_x):
    nodes = [i for i, value in enumerate(syndrom) if value == 1]
    graph = nx.complete_graph(nodes)
    for i,j in graph.edges:
        difference_rows = abs((i)//(n_qubits_x - 1) - (j)//(n_qubits_x - 1))
        max_val = max(i, j)
        min_val = min(i, j)
        graph[i][j]["weight"] = abs(max_val - difference_rows*(n_qubits_x - 1) - min_val) + difference_rows
    return graph

def visualise_graph(graph):
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos, with_labels=True)
    weights = nx.get_edge_attributes(graph, "weight")
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=weights)
    plt.show()

def create_correction(graph, n_qubits):
    list_correction = np.full(n_qubits, 'I')
    matching = nx.min_weight_matching(graph, weight="weight")
    for j, i in matching:
        list_correction

n_qubits_x = 7
n_qubits_y = 4
n_ancilla = n_qubits_x * n_qubits_y
qc = QuantumCircuit((n_qubits_x*(n_qubits_y + 1) + n_qubits_y**2) + n_ancilla, n_ancilla)
simulator = AerSimulator()
list_noise = create_noise(0.1 , None, n_qubits_x, n_qubits_y)
apply_noise(list_noise, qc)
syndrom = detection_syndrom(n_ancilla, qc, n_qubits_x, n_qubits_y, simulator)
qc.draw()
# syndrom = [0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0]
# graph = make_graph(syndrom , nx.Graph(), n_qubits_y, n_qubits_x)
# print(syndrom)
# visualise_graph(graph)

def visualise_lattice(n_qubits_x, n_qubits_y, syndrom, errors=None):
    n_data = n_qubits_x*(n_qubits_y + 1) + (n_qubits_x - 1)*n_qubits_y
    pos = {}   # indice du qubit de donnée -> (x, y)

    for i in range(n_qubits_x):
        # ligne du haut : n_y + 1 qubits
        for j in range(n_qubits_y + 1):
            pos[i*(2*n_qubits_y + 1) + j] = (2*j, 2*i)
        # ligne du dessous : n_y qubits
        if i < n_qubits_x - 1:
            for j in range(n_qubits_y):
                pos[i*(2*n_qubits_y + 1) + n_qubits_y + 1 + j] = (2*j + 1, 2*i + 1)

    fig, ax = plt.subplots(figsize=(9, 6))

    for q, (x, y) in pos.items():
        ax.plot(x, y, 'o', color='black', markersize=6, zorder=2)

    if errors is not None:
        for q in range(n_data):
            if errors[q] == 1:
                x, y = pos[q]
                ax.plot(x, y, 'o', color='slateblue', markersize=9, zorder=3)

    for i in range(n_qubits_x):
        for j in range(n_qubits_y):
            reference = i*(2*n_qubits_y + 1) + j
            x, y = 2*j + 1, 2*i

            voisins = [reference, reference + 1]
            if i > 0:
                voisins.append(reference - n_qubits_y)
            if i < n_qubits_x - 1:
                voisins.append(reference + n_qubits_y + 1)
            for v in voisins:
                vx, vy = pos[v]
                ax.plot([x, vx], [y, vy], color='gray', zorder=0)

            if syndrom[i*n_qubits_y + j] == 1:
                ax.plot(x, y, 'o', color='deepskyblue', markersize=22, zorder=4)
            else:
                ax.plot(x, y, 'o', color='lightgray', markersize=5, zorder=1)

    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.axis('off')
    plt.show()
visualise_lattice(n_qubits_x, n_qubits_y, syndrom)


