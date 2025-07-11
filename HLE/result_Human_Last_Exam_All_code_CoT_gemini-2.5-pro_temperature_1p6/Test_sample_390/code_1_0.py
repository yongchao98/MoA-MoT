import numpy as np

def solve():
    """
    This function demonstrates that for a set of orthogonal vectors {y_i},
    the set S forms a simplex.
    """
    # Define a set of 3 orthogonal vectors in R^4
    # They are linearly independent by definition.
    y1 = np.array([2, 0, 0, 0])
    y2 = np.array([0, 3, 0, 0])
    y3 = np.array([0, 0, 1, 0])
    
    y_vectors = [y1, y2, y3]
    n = len(y_vectors)
    
    # Calculate the norms squared ||y_i||^2
    y_norms_sq = [np.linalg.norm(y)**2 for y in y_vectors]
    
    # The equation of the simplex is sum(x_i / ||y_i||^2) = 1
    print("The vectors y_i are:")
    for i, y in enumerate(y_vectors):
        print(f"y{i+1} = {y}")
    
    print("\nThe squared norms ||y_i||^2 are:")
    for i, norm_sq in enumerate(y_norms_sq):
        print(f"||y{i+1}||^2 = {norm_sq}")
    
    print("\nThe equation of the simplex S is sum(w_i * x_i) = 1, where w_i = 1/||y_i||^2.")
    print("This specific equation is:")
    equation_parts = []
    for i in range(n):
        equation_parts.append(f"x_{i+1}/{y_norms_sq[i]:.2f}")
    print(" + ".join(equation_parts) + " = 1")
    
    print("\nLet's test this with 5 random points.")
    
    # Orthonormal basis for the span of {y_i}
    q, r = np.linalg.qr(np.array(y_vectors).T)
    basis = q[:,:n]

    for i in range(5):
        # Generate a random unit vector s in the span of {y_i}
        # A random vector in the span is a linear combination with random coefficients
        # of the basis vectors. Normalizing it makes it a unit vector.
        random_coeffs = np.random.randn(n)
        s = basis @ random_coeffs
        s = s / np.linalg.norm(s)
        
        # Calculate the point in S
        x = np.array([np.abs(np.dot(y, s))**2 for y in y_vectors])
        
        # Verify if the point lies on the simplex plane
        simplex_sum = np.sum(x / np.array(y_norms_sq))
        
        print(f"\n--- Test Point {i+1} ---")
        print(f"Random unit vector s (first 4 components): {s[:4]}")
        print(f"Corresponding point x in S: {x}")
        print(f"Check sum: {' + '.join([f'{xi:.4f}/{norm_sq:.2f}' for xi, norm_sq in zip(x, y_norms_sq)])} = {simplex_sum:.4f}")
        assert np.isclose(simplex_sum, 1.0), "Point does not lie on the simplex"

solve()