import math

def solve_quantum_circuit():
    """
    This function calculates and prints the step-by-step evolution of a 3-qubit quantum state
    and determines the final probability of measuring the |100> state.
    """
    
    # Define constants for printing
    s2 = "1/sqrt(2)"
    
    print("Step-by-step derivation of the final quantum state:\n")

    # Step 0: Initial State
    psi_0 = "|000>"
    print(f"Step 0: The initial state is")
    print(f"|ψ₀⟩ = {psi_0}\n")

    # Step 1: Apply Hadamard gate to the first qubit
    print("Step 1: Apply a Hadamard gate to the first qubit (H ⊗ I ⊗ I).")
    print(f"H|0⟩ = {s2}(|0⟩ + |1⟩)")
    psi_1_ket = f"{s2}(|000⟩ + |100⟩)"
    print(f"|ψ₁⟩ = (H ⊗ I ⊗ I) |ψ₀⟩ = {s2}(|0⟩ + |1⟩) ⊗ |00⟩ = {psi_1_ket}\n")

    # Step 2: Apply CNOT gate (control=1, target=2)
    print("Step 2: Apply a CNOT gate with qubit 1 as control and qubit 2 as target.")
    print(f"CNOT acts on each component of |ψ₁⟩:")
    print(f"  CNOT₁₂ |000⟩ = |000⟩")
    print(f"  CNOT₁₂ |100⟩ = |110⟩")
    psi_2_ket = f"{s2}(|000⟩ + |110⟩)"
    print(f"|ψ₂⟩ = {psi_2_ket}\n")

    # Step 3: Apply Toffoli gate (controls=1,2, target=3)
    print("Step 3: Apply a Toffoli (CCNOT) gate with qubits 1 and 2 as controls and qubit 3 as target.")
    print(f"CCNOT acts on each component of |ψ₂⟩:")
    print(f"  CCNOT₁₂₃ |000⟩ = |000⟩ (Controls 00 do not activate the gate)")
    print(f"  CCNOT₁₂₃ |110⟩ = |111⟩ (Controls 11 activate the gate, flipping the target from 0 to 1)")
    psi_3_ket = f"{s2}(|000⟩ + |111⟩)"
    print(f"|ψ₃⟩ = {psi_3_ket}\n")
    
    # Step 4: Apply second Hadamard gate to the first qubit
    print("Step 4: Apply a second Hadamard gate to the first qubit.")
    print(f"The gate acts on each component of |ψ₃⟩:")
    print(f"  (H ⊗ I ⊗ I) |000⟩ = (H|0⟩) ⊗ |00⟩ = {s2}(|000⟩ + |100⟩)")
    print(f"  (H ⊗ I ⊗ I) |111⟩ = (H|1⟩) ⊗ |11⟩ = {s2}(|011⟩ - |111⟩)")
    print(f"|ψ₄⟩ = {s2} [ {s2}(|000⟩ + |100⟩) + {s2}(|011⟩ - |111⟩) ]")
    psi_4_ket = "1/2 * (|000⟩ + |100⟩ + |011⟩ - |111⟩)"
    print(f"|ψ₄⟩ = {psi_4_ket}\n")
    
    # Final Calculation
    print("Final Probability Calculation:\n")
    print(f"The final state is |ψ₄⟩ = 1/2|000⟩ + 1/2|100⟩ + 1/2|011⟩ - 1/2|111⟩.")
    print("The probability of an outcome is the squared magnitude of its amplitude.")
    
    amplitude_val = 0.5
    amplitude_str = "1/2"
    
    print(f"The amplitude of the state |100⟩ is {amplitude_str}.")
    
    probability = amplitude_val ** 2
    
    # The prompt requires printing each number in the final equation.
    print(f"The probability P(|100⟩) = |{amplitude_str}|² = ({amplitude_val})² = {probability}")

solve_quantum_circuit()