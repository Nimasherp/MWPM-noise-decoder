import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from tqdm import tqdm

##########################################################################
# This code is a first hand into correction code for noises like bit-flip
# It is in 1D so the use of graphs and MWPM is not needed 
# It is a demo and a test code for the 2D program
##########################################################################


# "noise()" gives a list of noise, bit-flip or phase-flip, based on their probabilities
# each to be implemented on the qubits
def create_noise(list_noise, px, pz, n_qubits):
    for i in range(n_qubits):
        if np.random.random() < px:
            list_noise[i] = 'X'
        # if np.random.random() < pz:
        #     list_noise[i] = 'Z' if list_noise[i] == 'I' else 'Y'
    return list_noise

# noise_X_bool() focus on the bit-flip noise and make the reste of the noise False
def noise_X_bool(list_noise, n_qubits):
    bool_X = np.full(n_qubits, False)
    for i, bruit in enumerate(list_noise):
        if(bruit == 'X' or bruit == 'Y'):
            bool_X[i] = True
    return bool_X

def apply_noise(list_noise_bool, qc):
    for i, noise in enumerate(list_noise_bool):
        if noise :
            qc.x(i)
    return qc

# detection() is a function that will illuminate the syndrom by comparing the qubits next to eachother
# So in a real-life experiment we only get this information
def detection_syndrom(n_ancilla, qc, simulator):
    for i in range(n_ancilla):
        qc.cx(i, n_ancilla + 1 + i)
        qc.cx(i + 1, n_ancilla + i + 1)
    
    for i in range(n_ancilla):
        qc.measure(n_ancilla + 1 + i, i)
    
    result = simulator.run(qc, shots=1).result()

    bitstring = list(result.get_counts().keys())[0]
    syndrom = [int(bit) for bit in bitstring[::-1]]

    return syndrom

# make_graph creates thev complete graph of all syndroms
def make_graph(syndroms, graph):
    nodes = [i for i, value in enumerate(syndroms) if value == 1]
    graph = nx.complete_graph(nodes)
    for i,j in graph.edges:
        graph[i][j]["weight"] = abs(j - i)
    return graph

def apply_correction(graph, n_qubits):
    matching = nx.min_weight_matching(graph, weight="weight")
    list_correction = np.zeros(n_qubits, dtype=bool)
    for j, i in matching:
        list_correction[i + 1 : j + 1] = True
    return list_correction

def reset_ancilla(n_ancilla, qc):
    for i in range(n_ancilla):
        qc.reset(n_ancilla + 1 + i)
    return qc


def one_experiment(n_qubits, px, simulator):

    n_ancilla = n_qubits - 1

    list_noise = create_noise(np.full(n_qubits, 'I'), px, None, n_qubits)

    list_noise_bool = noise_X_bool(list_noise, n_qubits)

    qc = QuantumCircuit(n_qubits + n_ancilla, n_ancilla)

    

    apply_noise(list_noise_bool, qc)

    syndrome = detection_syndrom(n_ancilla, qc, simulator)

    # MWPM
    graph = make_graph(syndrome, nx.Graph())

    correction = apply_correction(graph, n_qubits)

    # correction
    reset_ancilla(n_ancilla, qc)
    apply_noise(correction, qc)

    error_before = np.sum(list_noise_bool)

    error_after = np.sum(
        list_noise_bool != correction
    )

    return error_before / n_qubits, error_after / n_qubits
    


def main():
    px_values = np.linspace(0.01, 0.6, 20)
    n_qubits = 20


    errors_before = []
    errors_after = []

    simulator = AerSimulator()

    n_tests = 1000

    for px in px_values:

        errors_before_px = []
        errors_after_px = []

        for _ in tqdm(
        range(n_tests),
        desc=f"pX = {px:.2f}"
        ):

            error_before, error_after = one_experiment(n_qubits, px, simulator)

            errors_before_px.append(error_before)
            errors_after_px.append(error_after)

        errors_before.append(np.mean(errors_before_px))
        errors_after.append(np.mean(errors_after_px))

    visualise_errors(errors_before, errors_after, px_values, n_qubits)

def visualise_errors(errors_before, errors_after, px_values, n_qubits):
    plt.figure(figsize=(8, 5))

    plt.plot(px_values, errors_before, marker='o', label="Before correction")

    plt.plot(px_values, errors_after, marker='o', label="After correction")

    plt.xlabel("Probability of bit-flip $p_X$")
    plt.ylabel("Rate of errors")

    plt.title(f"MWPM correction code ---- {n_qubits} qubits")
    plt.grid()
    plt.legend()
    plt.show()

main()

def visualise_graph(graph):
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos, with_labels=True)
    weights = nx.get_edge_attributes(graph, "weight")
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=weights)
    plt.show()


