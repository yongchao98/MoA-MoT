import sympy

def solve_peg_game_classes():
    """
    Calculates the number of equivalence classes in the peg game.

    The method is based on finding invariants using modular arithmetic.
    We consider the parity (even/odd) of the number of pegs for positions (x,y)
    grouped by their coordinates modulo 3. This gives 9 categories of positions,
    which we can think of as a 3x3 grid. Let p_ij be the parity for the
    category (x=i mod 3, y=j mod 3).

    A move (horizontal or vertical) flips the parities of three categories that
    form a row or a column in this 3x3 grid.

    An invariant is a linear combination of these 9 parities (mod 2) that is
    unaffected by any move. The coefficients c_ij of such a linear combination must
    satisfy a system of linear equations: the sum of coefficients for each
    row and each column must be 0 (mod 2).

    The number of independent invariants is the dimension of the solution space
    of this system. The total number of equivalence classes is 2 raised to the
    power of this dimension.
    """

    # There are 9 variables, c_ij, for the coefficients of the invariant.
    # We can flatten them into a 9-element vector:
    # [c_00, c_01, c_02, c_10, c_11, c_12, c_20, c_21, c_22]
    num_variables = 9

    # The constraints that sums of coefficients for rows and columns are 0 (mod 2)
    # form a system of linear equations Ac = 0.
    # The rows of matrix A represent these constraint equations.
    A_list = [
        # Row sum constraints
        [1, 1, 1, 0, 0, 0, 0, 0, 0],  # c_00 + c_01 + c_02 = 0
        [0, 0, 0, 1, 1, 1, 0, 0, 0],  # c_10 + c_11 + c_12 = 0
        [0, 0, 0, 0, 0, 0, 1, 1, 1],  # c_20 + c_21 + c_22 = 0
        # Column sum constraints
        [1, 0, 0, 1, 0, 0, 1, 0, 0],  # c_00 + c_10 + c_20 = 0
        [0, 1, 0, 0, 1, 0, 0, 1, 0],  # c_01 + c_11 + c_21 = 0
        [0, 0, 1, 0, 0, 1, 0, 0, 1],  # c_02 + c_12 + c_22 = 0
    ]

    # We use sympy to work with matrices over the finite field GF(2).
    A_matrix = sympy.Matrix(A_list)
    A_GF2 = A_matrix.applyfunc(lambda x: sympy.GF(2)(x))

    # The rank of the matrix tells us the number of independent constraints.
    rank = A_GF2.rank()

    # The dimension of the solution space (null space) is the number of variables
    # minus the rank of the matrix. This is our number of independent invariants.
    dim_invariants = num_variables - rank

    # The number of equivalence classes is 2 to the power of the number of invariants.
    num_classes = 2**dim_invariants

    print("Step 1: Define state by parities of peg counts in a 3x3 grid based on coordinates mod 3.")
    print(f"This gives {num_variables} parity variables (p_00, ..., p_22).")
    print("\nStep 2: Find the number of independent invariants.")
    print("An invariant is a combination of parities unchanged by any move.")
    print("This leads to a system of linear equations over GF(2).")
    print("The number of independent invariants is the dimension of the null space of the constraint matrix.")
    print(f"The constraint matrix has {A_matrix.shape[0]} rows (equations) and {A_matrix.shape[1]} columns (variables).")
    
    print(f"\nStep 3: Calculate the rank of the constraint matrix.")
    print(f"The rank over GF(2) is: {rank}")
    
    print("\nStep 4: Calculate the number of independent invariants.")
    print(f"Dimension of null space = (Number of variables) - (Rank)")
    print(f"Number of invariants = {num_variables} - {rank} = {dim_invariants}")

    print("\nStep 5: Calculate the total number of equivalence classes.")
    print("The number of classes is 2 to the power of the number of invariants.")
    print(f"Total number of classes = 2 ** {dim_invariants} = {num_classes}")

solve_peg_game_classes()