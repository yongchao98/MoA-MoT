import math

def solve_class_d_variables():
    """
    Calculates the number of non-Grassman variables for the bosonic sector of the
    supersymmetric sigma-model with N replicas for symmetry class D.
    """
    # Number of replicas as specified in the problem
    N = 2

    # The target space for the bosonic sector is the symmetric space SO(2N)/U(N).
    # The number of variables is the dimension of this space.
    # Dimension(G/H) = Dimension(G) - Dimension(H)

    # For G = SO(2N), the group order is 2*N
    group_order_G = 2 * N
    # Dimension of SO(n) = n * (n - 1) / 2
    dim_G = group_order_G * (group_order_G - 1) // 2

    # For H = U(N), the group order is N
    group_order_H = N
    # Dimension of U(n) = n^2
    dim_H = group_order_H ** 2

    # Calculate the final result
    num_variables = dim_G - dim_H

    # Print the explanation and the step-by-step calculation
    print(f"The number of non-Grassman (bosonic) variables for a class D SUSY sigma-model is the dimension of the symmetric space SO(2N)/U(N).")
    print(f"For this problem, the number of replicas N = {N}.")
    print("-" * 50)
    print(f"Step 1: Determine the dimension of the group G = SO(2*N) = SO({group_order_G}).")
    print(f"dim(SO(n)) = n * (n - 1) / 2")
    print(f"dim(SO({group_order_G})) = {group_order_G} * ({group_order_G} - 1) / 2 = {dim_G}")
    print("-" * 50)
    print(f"Step 2: Determine the dimension of the subgroup H = U(N) = U({group_order_H}).")
    print(f"dim(U(n)) = n^2")
    print(f"dim(U({group_order_H})) = {group_order_H}^2 = {dim_H}")
    print("-" * 50)
    print(f"Step 3: Calculate the dimension of the quotient space SO({group_order_G})/U({group_order_H}).")
    print(f"The dimension is dim(SO({group_order_G})) - dim(U({group_order_H})).")
    print(f"\nFinal Equation: {dim_G} - {dim_H} = {num_variables}")


if __name__ == '__main__':
    solve_class_d_variables()
