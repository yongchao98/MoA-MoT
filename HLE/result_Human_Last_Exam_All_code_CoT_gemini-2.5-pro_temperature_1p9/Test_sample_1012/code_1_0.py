# The accepted classical security level for long-term data protection is 128 bits.
# This means an attacker needs to perform 2^128 operations to break it.
classical_security_level_bits = 128

# Grover's algorithm, a quantum search algorithm, provides a quadratic speedup.
# This means the exponent of the work factor is effectively halved.
# To counteract this, the key length must be doubled to maintain the same security level.
quantum_resistance_factor = 2

# Calculate the minimum required key length for post-quantum security.
required_post_quantum_key_length_bits = classical_security_level_bits * quantum_resistance_factor

print("Step 1: Define the target classical security level.")
print(f"Target Security Level: {classical_security_level_bits} bits (requires 2^{classical_security_level_bits} operations to break).")
print("\nStep 2: Account for the quantum threat (Grover's Algorithm).")
print("Grover's algorithm effectively halves the bit-security of a key.")
print("To maintain the target security, we must double the key length.")
print("\nStep 3: Calculate the required post-quantum key length.")
print(f"Equation: Required Length = Target Security Level * Quantum Factor")
print(f"Calculation: {classical_security_level_bits} * {quantum_resistance_factor} = {required_post_quantum_key_length_bits} bits")

print(f"\nFinal Answer: A symmetric key needs to be at least {required_post_quantum_key_length_bits} bits long to be secure against a quantum computer.")
