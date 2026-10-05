import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from tqdm import tqdm
import itertools
from matplotlib.lines import Line2D
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
def create_correction(matching, graph_circuit, n_qubits, qc, n_ancilla, n_qubits_x):
    n_data = qc.num_qubits - n_ancilla
    n_x = n_qubits_x - 1
    list_correction = np.full(n_qubits, 'I')

    def real_border(syn_idx):
        # closest border
        row, col = divmod(syn_idx, n_x)
        side = 0 if col + 1 <= n_x - col else n_x - 1
        return ('B', n_data + row * n_x + side)

    for a, b in matching:
        a_is_border = isinstance(a, tuple)
        b_is_border = isinstance(b, tuple)

        if a_is_border and b_is_border:
            continue

        a = real_border(a[1]) if a_is_border else a + n_data
        b = real_border(b[1]) if b_is_border else b + n_data

        path = nx.shortest_path(graph_circuit, a, b)
        for node in path:
            if graph_circuit.nodes[node]["type"] == "qubit":
                list_correction[node] = 'X'

    return list_correction



def one_experiment(n_qubits_x, n_qubits_y, px, simulator):
    n_ancilla = (n_qubits_x - 1) * (n_qubits_y + 1)
    n_qubits = n_qubits_x*(n_qubits_y + 1) + n_qubits_y*(n_qubits_x - 1)

    qc = QuantumCircuit(n_qubits + n_ancilla, n_ancilla)
    
    list_noise = create_noise(px, None, n_qubits)
    error_before = np.count_nonzero(list_noise == 'X') / n_qubits
    apply_noise(list_noise, qc)
    n_data = qc.num_qubits - n_ancilla
    syndrom, graph_circuit = detection_syndrom(n_ancilla, qc, n_qubits_x, n_qubits_y,simulator, nx.Graph())

    
    matching = mwpm(syndrom, n_qubits_y, n_qubits_x)
    

    list_correction = create_correction(matching, graph_circuit, n_qubits, qc, n_ancilla, n_qubits_x)

    residual_error = np.where(
        list_noise == list_correction,
        'I',
        'X'
    )
    error_after = np.count_nonzero(residual_error == 'X') / n_qubits
    return error_before, error_after

def main():
    px_values = np.linspace(0.001, 0.1, 15)

    n_qubits_x = 25
    n_qubits_y = 20

    simulator = AerSimulator()

    n_tests = 20

    error_before = []
    error_after = []

    for px in px_values:

        before = 0
        after = 0

        for _ in tqdm(
        range(n_tests),
        desc=f"pX = {px:.2f}"
        ):

            result_before, result_after = one_experiment(
                n_qubits_x,
                n_qubits_y,
                px,
                simulator
            )

            before += result_before
            after += result_after

        error_before.append(before / n_tests)
        error_after.append(after / n_tests)


    plt.plot(px_values, error_before, "o-", label="Before MWPM")
    plt.plot(px_values, error_after, "o-", label="After MWPM")

    plt.xlabel("Bit-flip probability $p_x$")
    plt.ylabel("Fraction of qubits in error")
    plt.title("Physical error rate before and after MWPM")
    plt.legend()
    plt.grid(alpha=0.3)

    plt.show()
    
# one_experiment(4, 5, 0.4, AerSimulator())
main()



