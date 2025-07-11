import numpy as np

def analyze_shape_of_S(vectors):
    """
    Analyzes the shape of the set S for a given set of orthogonal vectors.

    Args:
        vectors (list of list of int/float): A list of vectors y_i.
                                              It's assumed these vectors are orthogonal.
    """
    try:
        # Convert to a numpy array for easier calculations
        y = np.array(vectors, dtype=float)
        n, d = y.shape
        
        # Verify that the vectors are linearly independent (rank of matrix is n)
        if np.linalg.matrix_rank(y) != n:
            print("Error: The provided vectors are not linearly independent.")
            return

        # Verify orthogonality
        gram_matrix = y @ y.T
        if not np.allclose(gram_matrix, np.diag(np.diag(gram_matrix))):
            print("Warning: The provided vectors are not orthogonal.")
            print("The analysis below assumes orthogonality and may not be correct for the general case.")

        # Calculate the squared norms (L_i)
        L = np.sum(y**2, axis=1)

        if np.any(L == 0):
            print("Error: One of the vectors is the zero vector, which is not allowed for linear independence.")
            return

        # The coefficients of the linear equation sum(c_i * x_i) = 1 are 1/L_i
        coeffs = 1.0 / L

        # Build and print the equation
        equation_parts = []
        for i in range(n):
            equation_parts.append(f"({coeffs[i]:.4f}) * x_{i+1}")

        equation = " + ".join(equation_parts) + " = 1"

        print("For the given orthogonal vectors, the set S forms a simplex.")
        print("The simplex lies on a hyperplane defined by the following linear equation:")
        print(equation)
        
        # Also print out each number in the equation individually as requested.
        print("\nThe numbers in the final equation are:")
        for coeff in coeffs:
            print(f"{coeff:.4f}")
        print("1")

    except Exception as e:
        print(f"An error occurred: {e}")

# Example with a set of 3 orthogonal vectors in R^3.
# y_1 = (2, 0, 0), y_2 = (0, 4, 0), y_3 = (0, 0, 1)
# L_1 = 4, L_2 = 16, L_3 = 1
# Equation should be (1/4)x_1 + (1/16)x_2 + (1/1)x_3 = 1
orthogonal_vectors = [
    [2, 0, 0],
    [0, 4, 0],
    [0, 0, 1]
]

analyze_shape_of_S(orthogonal_vectors)
