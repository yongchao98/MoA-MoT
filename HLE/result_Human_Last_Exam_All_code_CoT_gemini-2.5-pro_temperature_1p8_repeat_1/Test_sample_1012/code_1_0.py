def calculate_post_quantum_key_length():
    """
    Calculates the minimum symmetric key length needed to resist attacks
    from a quantum computer using Grover's algorithm.
    """
    # A widely accepted security level to protect against classical supercomputers is 128 bits.
    # An attacker would need to perform, on average, 2^128 operations.
    desired_security_level = 128

    # Grover's algorithm for quantum computers provides a quadratic speedup,
    # which means the number of operations is the square root of the classical case.
    # This effectively halves the bit strength of the key (e.g., security_bits = key_length / 2).
    # To find the required key length to achieve our desired security level, we must double it.
    grovers_factor = 2
    
    min_quantum_resistant_key_length = desired_security_level * grovers_factor

    print("### Post-Quantum Symmetric Key Length Calculation ###")
    print(f"Desired security level (against any computer): {desired_security_level} bits")
    print("Quantum computers can use Grover's algorithm, which reduces effective key strength by a factor of 2.")
    print("\nTo counteract this, the key length must be doubled.")
    print("The equation is: (New Key Length) / (Grover's Factor) = (Desired Security Level)")
    
    # Printing the final equation with the numbers filled in
    print(f"\nFinal Equation: {min_quantum_resistant_key_length} / {grovers_factor} = {desired_security_level}")

    print(f"\nConclusion: To achieve a security level of {desired_security_level} bits against a quantum computer, a symmetric key must be at least {min_quantum_resistant_key_length} bits long.")

calculate_post_quantum_key_length()