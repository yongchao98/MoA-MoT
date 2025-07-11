def calculate_quantum_safe_key_length():
    """
    Calculates the minimum symmetric key length required to be secure against
    a quantum computer attack using Grover's algorithm.
    """

    # The baseline security level we want to maintain, which is the current
    # standard for classical computers (e.g., AES-128).
    classical_security_level = 128

    # Grover's algorithm provides a quadratic speedup, effectively halving the
    # bit strength of a symmetric key. So, the key length must be doubled.
    quantum_impact_factor = 2

    # Calculate the required key length for quantum resistance.
    quantum_safe_key_length = classical_security_level * quantum_impact_factor

    print("To be secure against a quantum computer, a symmetric key's strength must be doubled.")
    print(f"The desired classical security level is {classical_security_level} bits.")
    print("The final calculation is:")
    print(f"{classical_security_level} * {quantum_impact_factor} = {quantum_safe_key_length}")
    print(f"\nTherefore, a symmetric key needs to be at least {quantum_safe_key_length} bits long.")

calculate_quantum_safe_key_length()