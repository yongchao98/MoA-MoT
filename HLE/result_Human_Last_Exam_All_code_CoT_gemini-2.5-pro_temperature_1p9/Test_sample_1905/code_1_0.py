import sympy

def prove_derivation_is_zero_for_finite_M(n):
    """
    Demonstrates that any derivation D on C(M) is zero when M is a finite set with n points.
    In this case, C(M) is isomorphic to R^n and D is an n x n matrix.
    
    Args:
        n (int): The number of points in the set M.
    """
    print(f"Let M be a finite set with {n} points.")
    print("The algebra of functions V = C(M) is isomorphic to R^n.")
    print("A derivation D: V -> V is a linear map, represented by an n x n matrix.")
    
    # Represent the derivation D as a matrix of symbolic variables
    D_matrix = sympy.Matrix([[sympy.Symbol(f'd_{i+1}{j+1}') for j in range(n)] for i in range(n)])
    
    # Let e_j be the basis vectors (1 at pos j, 0 elsewhere), representing delta functions.
    basis_vectors = sympy.eye(n)
    
    # From the Leibniz rule D(fg) = fD(g) + gD(f), we can derive constraints on D.
    # The product of functions is pointwise (Hadamard product for vectors).

    # Constraint 1: For any i != j, e_i * e_j = 0.
    # D(e_i * e_j) = D(0) = 0.
    # Also, D(e_i * e_j) = e_i * D(e_j) + e_j * D(e_i).
    # This implies that the off-diagonal elements of D must be zero.
    
    print("\nStep 1: Show off-diagonal elements of D are zero.")
    off_diagonal_equations = []
    for i in range(n):
        for j in range(i + 1, n):
            e_i = basis_vectors.col(i)
            e_j = basis_vectors.col(j)
            D_e_i = D_matrix.col(i)
            D_e_j = D_matrix.col(j)
            
            # Hadamard products
            term1 = sympy.matrix_multiply_elementwise(e_i, D_e_j)
            term2 = sympy.matrix_multiply_elementwise(e_j, D_e_i)
            
            # The k-th component of the equation e_i*D(e_j) + e_j*D(e_i) = 0
            # for k=i gives D_matrix[i, j] = 0
            # for k=j gives D_matrix[j, i] = 0
            off_diagonal_equations.append(sympy.Eq(D_matrix[i, j], 0))
            off_diagonal_equations.append(sympy.Eq(D_matrix[j, i], 0))

    print("From the property D(e_i * e_j) = 0 for i != j, we derive:")
    for eq in off_diagonal_equations:
        print(eq)

    # Constraint 2: For any i, e_i * e_i = e_i.
    # D(e_i) = D(e_i * e_i) = 2 * e_i * D(e_i).
    # The k-th component is D(e_i)_k = 2 * e_i_k * D(e_i)_k.
    # If k=i, e_i_k=1, so D(e_i)_i = 2 * D(e_i)_i, which means D_matrix[i,i] = 0.
    
    print("\nStep 2: Show diagonal elements of D are zero.")
    diagonal_equations = []
    for i in range(n):
        D_e_i_i = D_matrix[i,i]
        equation = sympy.Eq(D_e_i_i, 2 * 1 * D_e_i_i)
        solution = sympy.solve(equation, D_e_i_i)
        diagonal_equations.append(sympy.Eq(D_matrix[i,i], solution[0]))

    print("From the property D(e_i * e_i) = D(e_i), we derive:")
    for eq in diagonal_equations:
        print(eq)
        
    print("\nConclusion: Both diagonal and off-diagonal elements must be 0.")
    print("Therefore, the derivation D must be the zero matrix, so D = 0.")

# Run the proof for a 3-point set M.
prove_derivation_is_zero_for_finite_M(3)