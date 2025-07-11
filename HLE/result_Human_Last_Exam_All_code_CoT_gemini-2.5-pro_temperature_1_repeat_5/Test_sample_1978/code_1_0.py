import numpy as np

def solve_problem():
    """
    Calculates the index of the given boundary-value problem by determining the
    number of linearly independent boundary conditions.
    """
    # The index of a Boundary Value Problem is the number of linearly independent 
    # boundary conditions. We can represent the conditions as rows in a matrix 
    # and find its rank.
    # The variables involved are x_i(0) and x_i(T). We only need to consider
    # the variables that appear in the conditions: x1, x2, and x2024.
    # Let's create a coefficient matrix where columns correspond to:
    # x1(0), x2(0), x2024(0), x1(T), x2(T), x2024(T)

    # Original equations rewritten in the form a*x(0) + b*x(T) = 0:
    # 1: -1*x1(0) + 1*x1(T) = 0
    # 2: -1*x2(0) + 1*x2(T) = 0
    # 3: -5*x2(0) + 5*x2(T) = 0
    # 4: -100*x2(0) + 100*x2(T) = 0
    # 5: 1000*x2(0) - 1000*x2(T) = 0
    # 6: -100*x2024(0) + 100*x2024(T) = 0
    
    coefficient_matrix = np.array([
        [-1, 0, 0, 1, 0, 0],         # Equation 1
        [0, -1, 0, 0, 1, 0],         # Equation 2
        [0, -5, 0, 0, 5, 0],         # Equation 3 (5 * Equation 2)
        [0, -100, 0, 0, 100, 0],      # Equation 4 (100 * Equation 2)
        [0, 1000, 0, 0, -1000, 0],    # Equation 5 (-1000 * Equation 2)
        [0, 0, -100, 0, 0, 100]       # Equation 6
    ])

    # The rank of this matrix gives the number of linearly independent conditions.
    index = np.linalg.matrix_rank(coefficient_matrix)

    print("The boundary conditions can be simplified to a set of linearly independent equations.")
    print("The fundamental independent equations are:")
    # Equation 1
    print("1. For x1: {}*x1(0) + {}*x1(T) = 0".format(int(coefficient_matrix[0, 0]), int(coefficient_matrix[0, 3])))
    # Equation 2
    print("2. For x2: {}*x2(0) + {}*x2(T) = 0".format(int(coefficient_matrix[1, 1]), int(coefficient_matrix[1, 4])))
    # Equation 6
    print("3. For x2024: {}*x2024(0) + {}*x2024(T) = 0".format(int(coefficient_matrix[5, 2]), int(coefficient_matrix[5, 5])))
    
    print("\nThere are 3 linearly independent boundary conditions.")
    print(f"The index of the problem is the rank of the coefficient matrix, which is: {int(index)}")

solve_problem()