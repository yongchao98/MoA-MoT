def calculate_quantum_resistant_key_length():
    """
    Calculates the symmetric key length needed to resist a quantum computer attack.
    """
    # The current standard for strong security against classical computers is 128 bits.
    classical_security_level = 128

    print(f"To be secure against a classical computer, a symmetric key length of {classical_security_level} bits is considered a strong standard.")
    print("\nHowever, a quantum computer can use Grover's algorithm to speed up the key search.")
    print("This attack effectively halves the key's security strength in bits (n-bit key becomes n/2 bits of security).")
    print(f"To maintain a {classical_security_level}-bit security level against a quantum computer, the key length must be doubled.")

    # The quantum factor is 2 because Grover's algorithm offers a quadratic speedup.
    quantum_factor = 2
    quantum_resistant_key_length = classical_security_level * quantum_factor

    # Display the final equation with each number.
    print("\nThe calculation is as follows:")
    print(f"{classical_security_level} * {quantum_factor} = {quantum_resistant_key_length}")

    print(f"\nTherefore, a symmetric key needs to be at least {quantum_resistant_key_length} bits long to be considered secure against an arbitrarily powerful quantum computer.")

if __name__ == "__main__":
    calculate_quantum_resistant_key_length()