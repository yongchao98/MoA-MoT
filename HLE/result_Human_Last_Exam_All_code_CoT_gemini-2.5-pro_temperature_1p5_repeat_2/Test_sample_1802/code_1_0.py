import numpy as np

def calculate_max_rank_of_complementary_pauli_channel(d):
    """
    Calculates the maximal rank of the complementary channel of a d-dimensional
    Pauli channel by numerical verification.

    The maximal rank is achieved when all d^2 Pauli operators are included in
    the channel definition. The rank is the dimension of the space spanned
    by operator products U_ij^dagger * U_kl.

    Args:
        d (int): The dimension of the qudit system.

    Returns:
        int: The computed maximal rank.
    """
    if not isinstance(d, int) or d < 2:
        raise ValueError("Dimension 'd' must be an integer greater than or equal to 2.")

    # Define the shift (X) and clock (Z) matrices
    omega = np.exp(2j * np.pi / d)
    X = np.zeros((d, d), dtype=np.complex128)
    for i in range(d):
        X[i, (i - 1 + d) % d] = 1

    Z = np.diag([omega**i for i in range(d)])

    # Generate all d^2 generalized Pauli operators U_ij = X^i Z^j
    pauli_operators = []
    for i in range(d):
        for j in range(d):
            # Calculate U_ij = X^i * Z^j
            op = np.linalg.matrix_power(X, i) @ np.linalg.matrix_power(Z, j)
            pauli_operators.append(op)
    
    # Generate all d^4 products of the form U_k^dagger * U_l
    # and flatten them into vectors.
    product_vectors = []
    for uk in pauli_operators:
        for ul in pauli_operators:
            # Calculate the product and flatten it to a 1D vector
            product_matrix = uk.conj().T @ ul
            product_vectors.append(product_matrix.flatten())

    # Create a matrix where each row is a flattened operator product
    # The rank of this matrix is the dimension of the spanned space.
    basis_matrix = np.array(product_vectors)
    
    # The number of vectors can be large (d^4), but we only need d^2
    # linearly independent ones to span the space. Taking the rank of
    # this matrix composed of all products gives us the dimension of the span.
    rank = np.linalg.matrix_rank(basis_matrix)
    
    return int(rank)

# You can change the dimension 'd' here
d = 3

# Calculate and print the result
max_rank = calculate_max_rank_of_complementary_pauli_channel(d)
calculated_value = d**2

print(f"For a qudit of dimension d = {d}:")
print("The theoretical maximal rank of the complementary Pauli channel is d^2.")
print(f"The equation for the rank is: rank = {d}^2 = {calculated_value}")
print(f"The rank computed by the script is: {max_rank}")
if max_rank == calculated_value:
    print("\nThe computed value matches the theoretical result.")
else:
    print("\nThere is a discrepancy between the computed and theoretical result.")
