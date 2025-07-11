import numpy as np

def solve_pauli_channel_rank():
    """
    Calculates the maximal rank of the complementary channel of a d-dimensional
    Pauli channel by constructing the relevant vector space and finding its dimension.
    """
    # Let's choose a dimension 'd' for the qudit system as an example.
    # The result holds for any integer d >= 2. We'll use d=4.
    d = 4

    # The maximal rank of the complementary channel of a Pauli channel is d^2.
    # This is because the Choi matrix of a generic Pauli channel (with all non-zero
    # probabilities) is supported on a space spanned by d^2 orthogonal vectors.
    # These vectors are the generalized Bell basis states.
    
    # We will now numerically verify this by constructing these d^2 vectors
    # and computing the rank of the matrix formed by them.
    
    # Define the primitive d-th root of unity
    omega = np.exp(2j * np.pi / d)
    
    # The vectors that span the support of the Choi matrix are of the form:
    # |psi_ij> = sum_{k=0}^{d-1} omega^(j*k) |k> |(k+i) mod d>
    
    # Let's construct these d*d vectors. Each vector has d*d components.
    psi_vectors = []
    for i in range(d):
        for j in range(d):
            # Initialize a zero vector in the d*d dimensional space
            psi_ij = np.zeros(d*d, dtype=np.complex128)
            for k in range(d):
                # The basis state |k> |m> corresponds to index k*d + m
                # Here, m = (k+i) mod d
                index = k * d + (k + i) % d
                psi_ij[index] = omega**(j * k)
            psi_vectors.append(psi_ij)
            
    # Create a matrix M where the columns are the vectors |psi_ij>
    # The .T is because numpy treats list of lists as rows by default.
    M = np.array(psi_vectors).T

    # The number of linearly independent vectors is the rank of matrix M.
    maximal_rank = np.linalg.matrix_rank(M)
    
    # The final equation is Maximal Rank = d^2
    d_squared = d**2

    # Print out the results and the final equation
    print(f"For a qudit system of dimension d = {d}:")
    print(f"The maximal rank of the complementary channel is the dimension of the space spanned by the {d*d} generalized Bell basis vectors.")
    print(f"We construct these vectors and compute the rank of the matrix they form.")
    print("\n--- Calculation Result ---")
    print(f"The rank computed from the matrix is: {maximal_rank}")
    print(f"The theoretical value is d^2 = {d_squared}.")
    print("\n--- Final Equation ---")
    print(f"Maximal Rank = d * d = {d_squared}")

# Execute the function
solve_pauli_channel_rank()
