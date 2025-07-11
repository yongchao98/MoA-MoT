import math

def solve_quantum_circuit():
    """
    Solves the given quantum circuit problem step-by-step.
    """
    print("This script calculates the final state of a 3-qubit system and the probability of a specific measurement outcome.")
    print("-" * 80)

    # The state is represented as a dictionary mapping basis states (kets) to their amplitudes.
    # Initial State: |ψ₀⟩ = |000⟩
    state = {'000': 1.0}
    print("Step 0: Initial State")
    print("|ψ₀⟩ = |000⟩\n")

    # Step 1: Apply Hadamard gate to the first qubit (H ⊗ I ⊗ I)
    # H|0⟩ = 1/√2 (|0⟩ + |1⟩)
    # (H ⊗ I ⊗ I)|000⟩ = (H|0⟩) ⊗ |0⟩ ⊗ |0⟩ = 1/√2 (|0⟩ + |1⟩) ⊗ |00⟩ = 1/√2 (|000⟩ + |100⟩)
    psi_1_amp_000 = 1 / math.sqrt(2)
    psi_1_amp_100 = 1 / math.sqrt(2)
    state = {'000': psi_1_amp_000, '100': psi_1_amp_100}
    print("Step 1: Apply Hadamard gate to the first qubit.")
    print(f"|ψ₁⟩ = (H⊗I⊗I)|ψ₀⟩ = ({psi_1_amp_000:.3f})|000⟩ + ({psi_1_amp_100:.3f})|100⟩\n")

    # Step 2: Apply CNOT gate (control=1, target=2)
    # CNOT|000⟩ = |000⟩
    # CNOT|100⟩ = |110⟩ (flips 2nd qubit because 1st is 1)
    # |ψ₂⟩ = 1/√2 (CNOT|000⟩ + CNOT|100⟩) = 1/√2 (|000⟩ + |110⟩)
    psi_2_amp_000 = state['000']
    psi_2_amp_110 = state['100']
    state = {'000': psi_2_amp_000, '110': psi_2_amp_110}
    print("Step 2: Apply CNOT gate with qubit 1 as control and qubit 2 as target.")
    print(f"|ψ₂⟩ = CNOT₁₂|ψ₁⟩ = ({psi_2_amp_000:.3f})|000⟩ + ({psi_2_amp_110:.3f})|110⟩\n")

    # Step 3: Apply Toffoli gate (controls=1,2, target=3)
    # CCNOT|000⟩ = |000⟩ (controls 00, no flip)
    # CCNOT|110⟩ = |111⟩ (controls 11, flips 3rd qubit)
    # |ψ₃⟩ = 1/√2 (CCNOT|000⟩ + CCNOT|110⟩) = 1/√2 (|000⟩ + |111⟩)
    psi_3_amp_000 = state['000']
    psi_3_amp_111 = state['110']
    state = {'000': psi_3_amp_000, '111': psi_3_amp_111}
    print("Step 3: Apply Toffoli gate with qubits 1 and 2 as controls, and qubit 3 as target.")
    print(f"|ψ₃⟩ = CCNOT₁₂₃|ψ₂⟩ = ({psi_3_amp_000:.3f})|000⟩ + ({psi_3_amp_111:.3f})|111⟩\n")

    # Step 4: Apply a second Hadamard gate to the first qubit.
    # |ψ₄⟩ = (H⊗I⊗I)|ψ₃⟩ = 1/√2 * [ (H|0⟩)|00⟩ + (H|1⟩)|11⟩ ]
    # = 1/√2 * [ (1/√2)(|0⟩+|1⟩))|00⟩ + (1/√2)(|0⟩-|1⟩))|11⟩ ]
    # = 1/2 * [ |000⟩ + |100⟩ + |011⟩ - |111⟩ ]
    psi_4_amp_000 = state['000'] / math.sqrt(2)
    psi_4_amp_100 = state['000'] / math.sqrt(2)
    psi_4_amp_011 = state['111'] / math.sqrt(2)
    psi_4_amp_111 = -state['111'] / math.sqrt(2)
    state = {'000': psi_4_amp_000, '100': psi_4_amp_100, '011': psi_4_amp_011, '111': psi_4_amp_111}

    print("Step 4: Apply a second Hadamard gate to the first qubit.")
    print("The final state |ψ₄⟩ is composed of the following components:")
    print(f"|ψ₄⟩ = ({state['000']:.2f})|000⟩ + ({state['100']:.2f})|100⟩ + ({state['011']:.2f})|011⟩ + ({state['111']:.2f})|111⟩\n")

    # Final Step: Determine the probability of measuring |100⟩
    # The probability is the squared magnitude of the amplitude of |100⟩.
    amplitude_100 = state.get('100', 0)
    probability_100 = abs(amplitude_100)**2

    print("-" * 80)
    print("Final Calculation: Probability of measuring |100⟩")
    print(f"The amplitude of the |100⟩ component in |ψ₄⟩ is {amplitude_100:.2f}.")
    print("The probability is the squared magnitude of this amplitude.")
    print(f"P(|100⟩) = |amplitude|² = |{amplitude_100:.2f}|² = {probability_100:.4f}")

    print("\nThe final probability is 1/4.")
    print(f"Final answer as a decimal: {probability_100}")

solve_quantum_circuit()