import numpy as np

def solve_bvp_index():
    """
    Calculates the index of the given boundary-value problem.
    """
    # 1. Define the dimension of the system from x(t) = (x1, ..., x2024).
    n = 2024

    # 2. Define the number of given boundary conditions.
    m = 6

    # 3. Initialize the matrices B and C for the boundary conditions Bx(0) + Cx(T) = 0.
    # The matrices will have m rows (one for each condition) and n columns (one for each variable).
    B = np.zeros((m, n))
    C = np.zeros((m, n))

    # Populate B and C based on the given boundary conditions.
    # Note: Python uses 0-based indexing, so x_k corresponds to index k-1.

    # Condition 1: x_1(T) - x_1(0) = 0  => -1*x_1(0) + 1*x_1(T) = 0
    B[0, 0] = -1
    C[0, 0] = 1

    # Condition 2: x_2(T) - x_2(0) = 0  => -1*x_2(0) + 1*x_2(T) = 0
    B[1, 1] = -1
    C[1, 1] = 1

    # Condition 3: 5x_2(T) - 5x_2(0) = 0 => -5*x_2(0) + 5*x_2(T) = 0
    B[2, 1] = -5
    C[2, 1] = 5

    # Condition 4: 100x_2(T) - 100x_2(0) = 0 => -100*x_2(0) + 100*x_2(T) = 0
    B[3, 1] = -100
    C[3, 1] = 100

    # Condition 5: 1000x_2(0) - 1000x_2(T) = 0 => 1000*x_2(0) - 1000*x_2(T) = 0
    B[4, 1] = 1000
    C[4, 1] = -1000

    # Condition 6: 100x_2024(T) - 100x_2024(0) = 0 => -100*x_2024(0) + 100*x_2024(T) = 0
    B[5, 2023] = -100  # Index 2023 corresponds to x_2024
    C[5, 2023] = 100

    # 4. Form the characteristic matrix Q = B + C * e^(AT).
    # Since A is diagonal, e^(AT) is also diagonal. Let its diagonal elements be lambda_i.
    # The rank of Q is independent of the specific values of lambda_i, as long as they are not 1.
    # We can use a representative value, e.g., 0.5, for all lambda_i to compute the rank.
    lambdas = np.full(n, 0.5)
    exp_AT = np.diag(lambdas)
    
    # Calculate Q. We can do this efficiently without building the full diagonal matrix.
    Q = B + C * lambdas

    # 5. Calculate the rank of Q using numpy's linear algebra module.
    rank_Q = np.linalg.matrix_rank(Q)

    # 6. The index of the problem is n - rank(Q).
    index = n - int(rank_Q)

    # Print the details of the calculation.
    print(f"The dimension of the system is n = {n}.")
    print(f"The number of linearly independent boundary conditions is the rank of the characteristic matrix Q.")
    print(f"The rank of Q is {int(rank_Q)}.")
    print(f"The index is calculated by the formula: index = n - rank(Q).")
    print(f"index = {n} - {int(rank_Q)} = {index}")

if __name__ == "__main__":
    solve_bvp_index()