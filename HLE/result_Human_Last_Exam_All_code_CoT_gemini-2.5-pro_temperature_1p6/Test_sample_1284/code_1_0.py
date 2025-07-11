def find_smallest_dimension():
    """
    This function finds the smallest dimension 'n' for which the given Fourier extension
    inequality is known to fail.

    The explanation for the failure rests on a bilinear counterexample, which requires the
    dimension of the parameter space of the paraboloid (n-1) to be at least 2.
    """

    # The minimal dimension of the parameter space for the counterexample to work.
    min_param_dim = 2

    # The dimension of the parameter space is n-1. So, we solve for the smallest integer n in:
    # n - 1 >= min_param_dim
    # The smallest integer n satisfying this is n = min_param_dim + 1.
    smallest_n = min_param_dim + 1

    print("The condition for the counterexample to exist is that the dimension of the parameter space, n-1, must be at least 2.")
    print(f"Let n_minus_1 = {min_param_dim}")
    print(f"The smallest integer dimension n is n_minus_1 + 1 = {smallest_n}")


find_smallest_dimension()
