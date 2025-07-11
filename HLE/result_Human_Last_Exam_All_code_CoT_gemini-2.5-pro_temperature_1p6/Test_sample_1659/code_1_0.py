import math

def calculate_quantum_probability():
    """
    This function follows the evolution of a 3-qubit state through a quantum circuit
    and calculates the probability of measuring a specific outcome.
    """

    print("The evolution of the quantum state is as follows:")
    print("Initial state: |ψ₀> = |000>\n")

    # Step 1: Apply Hadamard gate to the first qubit
    # H|0> = 1/sqrt(2) * (|0> + |1>)
    # So, |ψ₁> = (H ⊗ I ⊗ I) |000> = 1/sqrt(2) * (|000> + |100>)
    print("1. After applying Hadamard gate to qubit 1:")
    print("|ψ₁> = (1/√2)|000> + (1/√2)|100>\n")

    # Step 2: Apply CNOT gate with qubit 1 as control, qubit 2 as target
    # CNOT|000> = |000>
    # CNOT|100> = |110>
    # So, |ψ₂> = 1/sqrt(2) * (|000> + |110>)
    print("2. After applying CNOT gate (control=1, target=2):")
    print("|ψ₂> = (1/√2)|000> + (1/√2)|110>\n")

    # Step 3: Apply Toffoli gate with qubits 1,2 as control, qubit 3 as target
    # CCNOT|000> = |000> (controls 00 do not flip target)
    # CCNOT|110> = |111> (controls 11 flip target)
    # So, |ψ₃> = 1/sqrt(2) * (|000> + |111>)
    print("3. After applying Toffoli gate (controls=1,2, target=3):")
    print("|ψ₃> = (1/√2)|000> + (1/√2)|111>\n")

    # Step 4: Apply Hadamard gate to the first qubit again
    # H|0> = 1/√2 * (|0> + |1>) -> H|000> = 1/√2 * (|000> + |100>)
    # H|1> = 1/√2 * (|0> - |1>) -> H|111> = 1/√2 * (|011> - |111>)
    # So, |ψ₄> = 1/√2 * ( 1/√2 * (|000> + |100>) + 1/√2 * (|011> - |111>) )
    # |ψ₄> = 1/2 * (|000> + |100> + |011> - |111>)
    print("4. After applying the final Hadamard gate to qubit 1:")

    # Define the final state's amplitudes
    c000 = 0.5
    c100 = 0.5
    c011 = 0.5
    c111 = -0.5
    
    print("The final state equation is:")
    print(f"|ψ₄> = ({c000})|000> + ({c100})|100> + ({c011})|011> + ({c111})|111>\n")

    # Probability Calculation
    amplitude_to_measure = c100
    probability = abs(amplitude_to_measure)**2

    print("The probability of measuring the outcome |100> is the squared magnitude of its amplitude.")
    print(f"The amplitude of |100> in the final state |ψ₄> is {amplitude_to_measure}.")
    print(f"Probability(|100>) = |{amplitude_to_measure}|² = {probability}")

if __name__ == '__main__':
    calculate_quantum_probability()