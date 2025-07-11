import sympy

def check_continuity_points(matrix_list):
    """
    Checks if matrices in a list are points of continuity for the minimal
    polynomial map.

    A matrix M is a point of continuity if and only if it is non-derogatory,
    which is equivalent to the degree of its minimal polynomial being equal to n.

    Args:
        matrix_list: A list of matrices, where each matrix is a list of lists.
    """
    print("--- Checking matrices for continuity property ---")
    for i, M_raw in enumerate(matrix_list):
        try:
            M = sympy.Matrix(M_raw)
            n = M.shape[0]
            if n == 0:
                print(f"\nMatrix {i+1}: [] is an empty matrix.")
                continue

            # The variable for the polynomial, can be any symbol.
            x = sympy.symbols('x')
            
            # Compute the minimal polynomial
            min_poly = M.minpoly(x)
            
            # Get the degree of the minimal polynomial
            deg_min_poly = sympy.degree(min_poly, gen=x)
            
            print(f"\nMatrix {i+1}:")
            sympy.pprint(M)
            print(f"\nMinimal polynomial: {min_poly}")
            # Final equation part from prompt: Explicitly show the comparison
            print(f"Degree of minimal polynomial: {deg_min_poly}")
            print(f"Size of matrix (n): {n}")
            
            if deg_min_poly == n:
                print(f"Result: Degree ({deg_min_poly}) == n ({n}). This is a point of continuity.")
            else:
                print(f"Result: Degree ({deg_min_poly}) != n ({n}). This is NOT a point of continuity.")
        except Exception as e:
            print(f"\nCould not process matrix {i+1}: {M_raw}. Error: {e}")
    print("\n--- Check complete ---")

if __name__ == '__main__':
    # Example matrices for testing
    
    # M1: 2x2 matrix with distinct eigenvalues (1, 2). Non-derogatory.
    M1 = [[0, 1], [-2, 3]]
    
    # M2: 2x2 scalar matrix. Derogatory. Minimal poly is x-1, degree 1 < 2.
    M2 = [[1, 0], [0, 1]]
    
    # M3: 3x3 matrix. Non-derogatory. Jordan form is J_2(1) + J_1(2).
    # Minimal polynomial is (x-1)^2 * (x-2), degree 3.
    M3 = [[1, 1, 0], [0, 1, 0], [0, 0, 2]]
    
    # M4: 3x3 matrix. Derogatory. Jordan form is J_1(1) + J_1(1) + J_1(2).
    # Minimal polynomial is (x-1)*(x-2), degree 2 < 3.
    M4 = [[1, 0, 0], [0, 1, 0], [0, 0, 2]]

    matrices_to_test = [M1, M2, M3, M4]
    check_continuity_points(matrices_to_test)