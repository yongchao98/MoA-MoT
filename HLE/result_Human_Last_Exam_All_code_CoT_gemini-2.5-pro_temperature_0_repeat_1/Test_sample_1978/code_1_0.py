import numpy as np

def solve_problem_index():
    """
    This function calculates the index of the given boundary-value problem
    by determining the number of linearly independent boundary conditions.
    """
    # The boundary conditions can be written in a matrix form where each row
    # represents the coefficients of a single condition. The variables involved
    # are x_1, x_2, and x_2024 at times t=0 and t=T.
    # We create a matrix where columns correspond to the coefficients of:
    # [x_1(0), x_1(T), x_2(0), x_2(T), x_2024(0), x_2024(T)]

    # The given boundary conditions are:
    # 1. x_1(T) - x_1(0) = 0              => -1*x_1(0) + 1*x_1(T) = 0
    # 2. x_2(T) - x_2(0) = 0              => -1*x_2(0) + 1*x_2(T) = 0
    # 3. 5*x_2(T) - 5*x_2(0) = 0          => -5*x_2(0) + 5*x_2(T) = 0
    # 4. 100*x_2(T) - 100*x_2(0) = 0      => -100*x_2(0) + 100*x_2(T) = 0
    # 5. 1000*x_2(0) - 1000*x_2(T) = 0    => 1000*x_2(0) - 1000*x_2(T) = 0
    # 6. 100*x_2024(T) - 100*x_2024(0) = 0 => -100*x_2024(0) + 100*x_2024(T) = 0

    coefficient_matrix = np.array([
        [-1, 1, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0],
        [0, 0, -5, 5, 0, 0],
        [0, 0, -100, 100, 0, 0],
        [0, 0, 1000, -1000, 0, 0],
        [0, 0, 0, 0, -100, 100]
    ])

    # The index of the problem is the rank of this coefficient matrix.
    problem_index = np.linalg.matrix_rank(coefficient_matrix)

    # From our analysis, we have:
    # - 1 independent condition for x_1.
    # - 1 independent condition for x_2 (as conditions 2, 3, 4, 5 are dependent).
    # - 1 independent condition for x_2024.
    num_cond_x1 = 1
    num_cond_x2 = 1
    num_cond_x2024 = 1
    
    print("The index of the problem is the total number of linearly independent boundary conditions.")
    print(f"Number of independent conditions for x_1: {num_cond_x1}")
    print(f"Number of independent conditions for x_2: {num_cond_x2}")
    print(f"Number of independent conditions for x_2024: {num_cond_x2024}")
    
    total_index = num_cond_x1 + num_cond_x2 + num_cond_x2024
    
    print("\nThe final equation for the index is:")
    # Printing the equation as requested
    print(f"{num_cond_x1} + {num_cond_x2} + {num_cond_x2024} = {total_index}")

    print("\nThe numbers in the final equation are:")
    print(num_cond_x1)
    print(num_cond_x2)
    print(num_cond_x2024)
    print(total_index)

    print(f"\nThis result is confirmed by calculating the rank of the coefficient matrix, which is: {problem_index}")

solve_problem_index()