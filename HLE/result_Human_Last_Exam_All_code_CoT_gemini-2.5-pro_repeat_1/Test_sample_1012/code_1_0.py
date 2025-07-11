# The current standard for strong classical security is 128 bits.
classical_security_level = 128

# Grover's algorithm provides a quadratic speedup, effectively halving the bit-strength of a symmetric key.
# To maintain the desired security level against a quantum computer, we must double the key length.
quantum_multiplier = 2

# Calculate the required key length for post-quantum security.
post_quantum_key_length = classical_security_level * quantum_multiplier

print("To achieve a security level of 128 bits against a quantum computer, a symmetric key must be at least:")
print(f"{classical_security_level} * {quantum_multiplier} = {post_quantum_key_length} bits")
