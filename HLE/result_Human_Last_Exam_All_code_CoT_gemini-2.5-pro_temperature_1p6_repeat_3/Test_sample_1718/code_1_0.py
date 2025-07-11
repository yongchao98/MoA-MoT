def solve_task():
    """
    Calculates the number of independent components of the Riemann tensor
    for a Kähler manifold of a given complex dimension.
    """
    # The real dimension of a Kähler manifold is even, n = 2m.
    # We will use the complex dimension 'm' for the calculation.
    # Let's calculate for a non-trivial example, a manifold of complex dimension m=2
    # (i.e., a real 4-dimensional manifold).
    m = 2

    # The number of independent components is determined by the number of
    # independent components of a Hermitian form on the space of symmetric
    # (2,0)-tensors. The dimension of this space is N = m(m+1)/2.
    N_dim_numerator = m * (m + 1)
    N_dim = N_dim_numerator // 2

    # The number of independent real components of an N x N Hermitian form is N^2.
    result = N_dim ** 2

    # As requested, we print each number in the final equation.
    print(f"For a Kähler manifold of complex dimension m = {m}:")
    print("The number of independent components is given by (m * (m + 1) / 2)^2.")
    print("The calculation is:")
    # Print the equation step-by-step
    equation = f"({m} * ({m} + 1) / 2)^2 = ({m} * {m+1} / 2)^2 = ({N_dim_numerator} / 2)^2 = ({N_dim})^2 = {result}"
    print(equation)

solve_task()