def solve_quantum_symmetric_key_length():
    """
    Calculates the minimum symmetric key length to be secure against a quantum computer.
    """
    classical_security_level = 128

    print("To determine the required symmetric key length for quantum resistance, we need to account for Grover's algorithm.")
    print("-" * 70)
    print(f"1. A widely accepted standard for strong security against classical computers is a {classical_security_level}-bit key.")
    print(f"   This implies that performing 2^{classical_security_level} operations is considered computationally infeasible.")
    print("\n2. Grover's algorithm, a quantum search algorithm, provides a quadratic speedup for breaking symmetric keys.")
    print("   For a key of length 'n', a quantum computer needs sqrt(2^n) = 2^(n/2) operations.")
    print("\n3. To be 'quantum-secure', the effort for a quantum computer must be at least as difficult as our classical security baseline.")
    print("   Therefore, we set the quantum effort equal to the classical security target:")
    print("   Quantum Effort = Classical Security Target")
    print("   2^(n/2) = 2^128")
    print("\n4. By comparing the exponents, we can solve for the new required key length 'n':")
    print("   n / 2 = 128")

    # Perform the final calculation
    required_n = classical_security_level * 2

    print(f"   n = {classical_security_level} * 2")
    print(f"   n = {required_n}")
    print("-" * 70)
    print(f"\nConclusion: A symmetric key must be at least {required_n} bits long to provide a security level of {classical_security_level} bits against a quantum computer.")

solve_quantum_symmetric_key_length()