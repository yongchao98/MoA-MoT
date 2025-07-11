# The desired security level in bits. A 128-bit security level is the
# current industry standard, meaning it requires 2^128 operations to break.
desired_security_level = 128

# Grover's algorithm, a quantum search algorithm, provides a quadratic
# speedup. This means it can break a symmetric key of length 'k' in
# roughly sqrt(2^k) = 2^(k/2) operations.

# To maintain our desired security level against a quantum computer,
# the post-quantum effort must be at least 2^128.
# So, we set the quantum attack complexity equal to the desired security level:
# k / 2 = desired_security_level

required_key_length = desired_security_level * 2

print(f"To achieve a security level of {desired_security_level} bits against a quantum computer,")
print("we must account for the quadratic speedup from Grover's algorithm.")
print("The relationship is: Post-Quantum Key Length / 2 = Desired Security Level.")
print(f"Therefore, the calculation is: {required_key_length} / 2 = {desired_security_level}")
print(f"The minimum required symmetric key length is {required_key_length} bits.")