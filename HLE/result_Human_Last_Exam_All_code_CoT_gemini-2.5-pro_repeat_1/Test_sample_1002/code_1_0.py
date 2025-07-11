import sys

def calculate_limit(k):
    """
    Calculates the value of the limit lim_{m->inf} ln(f(m))/ln(m).

    Args:
        k (int): An integer greater than or equal to 2.

    Returns:
        float: The value of the limit.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.", file=sys.stderr)
        return None

    # The formula for the limit is 1 - 1/(2k)
    numerator = 2 * k - 1
    denominator = 2 * k
    result = numerator / denominator

    # Outputting the final equation as requested
    print(f"For k = {k}, the final equation for the limit is:")
    print(f"1 - 1 / (2 * {k}) = {numerator} / {denominator}")
    print(f"The result is: {result}")

    return result

if __name__ == '__main__':
    # You can change the value of k here.
    # The problem states k is an integer >= 2.
    k_value = 3
    calculate_limit(k_value)

    print("\n--- Another example ---")
    k_value_2 = 4
    calculate_limit(k_value_2)
