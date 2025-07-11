import math

def solve_class_d_variables():
    """
    Calculates the number of non-Grassman variables for the bosonic sector of the
    supersymmetric sigma-model for class D with a given number of replicas.
    """
    # Number of replicas as specified in the problem
    n = 2

    # The number of non-Grassman (bosonic) variables corresponds to the dimension
    # of the bosonic part of the target supermanifold for the sigma-model.
    # For symmetry class D, this manifold is a product of two symmetric spaces:
    # M_B = [O(2n) / (O(n) x O(n))] x [Sp(4n) / (Sp(2n) x Sp(2n))]
    # where 'n' is the number of replicas.

    # We calculate the dimension of each part.

    # Part 1: Non-compact sector O(2n) / (O(n) x O(n))
    # The dimension of this space is n^2.
    dim_non_compact = n**2

    # Part 2: Compact sector Sp(4n) / (Sp(2n) x Sp(2n))
    # The dimension of this space is 4*n^2.
    dim_compact = 4 * n**2

    # The total number of variables is the sum of the dimensions.
    total_dim = dim_non_compact + dim_compact

    # Print the explanation and the final calculation
    print("The supersymmetric sigma-model for disordered systems of symmetry class D with n replicas")
    print("has a bosonic sector parametrized by a manifold with a dimension calculated as follows:")
    print(f"For n = {n} replicas:")
    print(f"Dimension of non-compact part (O({2*n})/(O({n})xO({n}))): n^2 = {n}^2 = {dim_non_compact}")
    print(f"Dimension of compact part (Sp({4*n})/(Sp({2*n})xSp({2*n}))): 4*n^2 = 4*{n}^2 = {dim_compact}")
    print("\nThe total number of non-Grassman (bosonic) variables is the sum of these dimensions.")
    print(f"Total variables = {dim_non_compact} + {dim_compact} = {total_dim}")

solve_class_d_variables()