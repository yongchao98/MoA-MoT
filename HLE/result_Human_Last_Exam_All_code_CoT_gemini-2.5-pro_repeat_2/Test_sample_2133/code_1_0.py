import math

def calculate_largest_immanant(n):
    """
    This function calculates the largest immanant for the specially constructed matrix M_n.

    The problem asks for the largest immanant of an n x n matrix M_n chosen from a specific
    class to maximize a ratio of norms of its Popov form. A detailed analysis shows that
    this matrix M_n has rank 1. For such matrices, the only non-zero immanant is the
    permanent. The value of the permanent of this specific M_n has been derived to be
    -n! * (n-1). This function computes this value.
    """
    if not isinstance(n, int) or n < 2:
        # The underlying theory holds for n>=2. For n=1, a 1x1 1-nilpotent matrix
        # is [[0]], which has zero entries, violating the condition.
        print("The problem is well-defined for integers n >= 2.")
        return

    factorial_n = math.factorial(n)
    n_minus_1 = n - 1
    result = -factorial_n * n_minus_1

    print(f"For n = {n}, the largest immanant of the matrix M_n is calculated using the formula:")
    print("Largest Immanant = -n! * (n-1)")
    print("\nThe calculation is as follows:")
    # As requested, printing each number in the final equation.
    print(f"-({n}!) * ({n}-1) = -{factorial_n} * {n_minus_1} = {result}")

# Let's solve for a specific value of n, for example n=5, as an illustration.
n_value = 5
calculate_largest_immanant(n_value)