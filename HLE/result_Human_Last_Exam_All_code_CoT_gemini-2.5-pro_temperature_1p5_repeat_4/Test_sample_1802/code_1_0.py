import numpy as np

def calculate_maximal_rank_of_complementary_pauli_channel(d):
    """
    This function calculates the maximal rank of the complementary channel of a
    d-dimensional Pauli channel.

    The rank is shown to be d**2. This code verifies it numerically by
    constructing the d**2 generalized Pauli operators and checking that they
    are linearly independent by computing the rank of the matrix formed
    by their vectorized representations.

    Args:
        d (int): The dimension of the qudit system.
    """
    if not isinstance(d, int) or d < 2:
        print("Error: Dimension 'd' must be an integer greater than or equal to 2.")
        return

    # Generate the generalized Pauli X (shift) and Z (phase) operators
    X = np.zeros((d, d), dtype=complex)
    for i in range(d):
        X[i, (i - 1 + d) % d] = 1

    omega = np.exp(2 * np.pi * 1j / d)
    Z = np.diag([omega**i for i in range(d)])

    # The maximal rank is achieved when all Pauli operators are in the basis
    # of the output space of the complementary channel.
    # We check the dimension of the space spanned by these operators.
    # We generate all d**2 Pauli operators U_ab = X^a Z^b
    pauli_operators = []
    for a in range(d):
        for b in range(d):
            # U_ab = (X^a) @ (Z^b)
            # It's better to compute as Z^b X^a with a phase, but for linear
            # independence, the order or phase doesn't matter.
            # We use U_ab = X^a Z^b for simplicity.
            Xa = np.linalg.matrix_power(X, a)
            Zb = np.linalg.matrix_power(Z, b)
            U_ab = Xa @ Zb
            pauli_operators.append(U_ab)

    # To find the dimension of the span, we vectorize each d x d operator
    # into a (d*d)-dimensional vector and form a matrix where these
    # vectors are the columns.
    vectorized_operators = []
    for op in pauli_operators:
        vectorized_operators.append(op.flatten())

    # Create a matrix from the list of vectors.
    # The matrix has d**2 rows and d**2 columns.
    matrix_to_test = np.array(vectorized_operators).T

    # The rank of this matrix gives the number of linearly independent operators,
    # which is the dimension of the space they span.
    rank = np.linalg.matrix_rank(matrix_to_test)
    
    # The analytical result for the maximal rank is d^2.
    maximal_rank = d**2
    
    # The final equation is: maximal_rank = d^2
    # We output each number and symbol in this equation.
    print(f"For a qudit of dimension d = {d}:")
    print("The theoretical maximal rank of the complementary Pauli channel is d^2.")
    print(f"Equation: maximal_rank = {d}^{2} = {maximal_rank}")
    print(f"Numerical verification of the rank: {rank}")
    if rank == maximal_rank:
        print("The numerical result matches the theoretical prediction.")
    else:
        print("The numerical result does NOT match the theoretical prediction.")

if __name__ == '__main__':
    # You can change this value to test for different dimensions.
    dimension = 3
    calculate_maximal_rank_of_complementary_pauli_channel(dimension)
