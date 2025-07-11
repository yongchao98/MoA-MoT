import math

def solve_quantum_circuit():
    """
    This function calculates and explains the step-by-step evolution of a 3-qubit state
    through a quantum circuit and determines the final probability of measuring |100>.
    """
    
    print("This script will determine the final state of a 3-qubit system and the probability of a specific measurement outcome.")
    print("-" * 80)
    
    # Initial State
    print("Initial State: The system starts in the state |ψ₀> = |000>.\n")
    
    # Step 1: H gate on qubit 1
    print("Step 1: Apply a Hadamard gate to the first qubit.")
    print("The operation is |ψ₁> = (H ⊗ I ⊗ I) |ψ₀> = H|0> ⊗ I|0> ⊗ I|0>.")
    print(f"Since H|0> = (1/√2)(|0> + |1>) = ({1/math.sqrt(2):.4f})(|0> + |1>), the state becomes:")
    print("|ψ₁> = (1/√2)(|0> + |1>) ⊗ |0> ⊗ |0>")
    print("|ψ₁> = (1/√2)(|000> + |100>)\n")
    
    # Step 2: CNOT gate
    print("Step 2: Apply a CNOT gate with qubit 1 as control and qubit 2 as target.")
    print("The operation is |ψ₂> = CNOT₁,₂ |ψ₁> = CNOT₁,₂ [(1/√2)(|000> + |100>)].")
    print("|ψ₂> = (1/√2) * [CNOT₁,₂|000> + CNOT₁,₂|100>]")
    print("The CNOT₁,₂ gate flips the second qubit if the first is |1>. So, CNOT₁,₂|000> = |000> and CNOT₁,₂|100> = |110>.")
    print("|ψ₂> = (1/√2)(|000> + |110>)\n")
    
    # Step 3: Toffoli gate
    print("Step 3: Apply a Toffoli (CCNOT) gate with qubits 1 and 2 as controls and qubit 3 as target.")
    print("The operation is |ψ₃> = CCNOT₁,₂,₃ |ψ₂> = CCNOT₁,₂,₃ [(1/√2)(|000> + |110>)].")
    print("|ψ₃> = (1/√2) * [CCNOT₁,₂,₃|000> + CCNOT₁,₂,₃|110>]")
    print("The CCNOT gate flips the third qubit if the first two are |1>. So, CCNOT|000> = |000> and CCNOT|110> = |111>.")
    print("|ψ₃> = (1/√2)(|000> + |111>)\n")
    
    # Step 4: H gate on qubit 1
    print("Step 4: Apply a second Hadamard gate to the first qubit.")
    print("The operation is |ψ₄> = (H ⊗ I ⊗ I) |ψ₃> = (H ⊗ I ⊗ I) [(1/√2)(|000> + |111>)].")
    print("|ψ₄> = (1/√2) * [(H⊗I⊗I)|000> + (H⊗I⊗I)|111>]")
    print("We have H|0> = (1/√2)(|0> + |1>) and H|1> = (1/√2)(|0> - |1>).")
    print("(H⊗I⊗I)|000> = (H|0>)⊗|00> = (1/√2)(|000> + |100>)")
    print("(H⊗I⊗I)|111> = (H|1>)⊗|11> = (1/√2)(|011> - |111>)")
    print("Substituting these back gives:")
    print("|ψ₄> = (1/√2) * [ (1/√2)(|000> + |100>) + (1/√2)(|011> - |111>) ]")
    print("|ψ₄> = 1/2 * (|000> + |100> + |011> - |111>)\n")

    # Final Calculation
    print("-" * 80)
    print("Calculation of the Probability")
    print("-" * 80)
    print("The final state is |ψ₄> = (1/2)|000> + (1/2)|100> + (1/2)|011> - (1/2)|111>.")
    print("The probability of measuring an outcome is the squared magnitude of its amplitude.")
    
    amplitude_str = "1/2"
    amplitude_val = 0.5
    probability = amplitude_val**2
    
    print(f"The amplitude of the outcome |100> in the final state |ψ₄> is {amplitude_str}.")
    print(f"Therefore, the probability P(|100>) is |{amplitude_str}|².")
    print(f"P(|100>) = |{amplitude_val}|² = {probability}")

solve_quantum_circuit()
<<<0.25>>>