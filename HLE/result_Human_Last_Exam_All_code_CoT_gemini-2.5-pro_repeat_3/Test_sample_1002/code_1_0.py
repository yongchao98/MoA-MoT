import sys

def solve_limit(k):
    """
    Calculates the value of the limit lim_{m->inf} (ln f(m) / ln m) for a given k.

    Args:
        k (int): An integer greater than or equal to 2.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: Please provide an integer k >= 2.", file=sys.stderr)
        return

    # The formula for the limit is 1 - 1/(2*k)
    numerator = 1
    denominator = 2 * k
    result = 1 - numerator / denominator

    print(f"For k = {k}, the problem is to find the limit of ln(f(m))/ln(m).")
    print(f"The derived formula for the limit is: 1 - 1 / (2 * k)")
    print(f"Plugging in k = {k}, the final expression is: 1 - {numerator} / {denominator}")
    print(f"The computed value of the limit is: {result}")


# You can change the value of k here.
# For example, let's use k=2 as a sample case.
k_value = 2
solve_limit(k_value)

# Example with k=3
# k_value = 3
# solve_limit(k_value)
