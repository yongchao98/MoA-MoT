import sys

def compute_limit(k):
    """
    This function computes the value of the limit based on the derived formula.
    The problem is to compute lim_{m -> infinity} (ln f(m)) / (ln m).
    Based on extremal graph theory, this limit evaluates to (k-1)/(k+1).

    Args:
        k (int): An integer k >= 2.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.", file=sys.stderr)
        return

    # The formula for the limit is (k-1)/(k+1)
    numerator = k - 1
    denominator = k + 1
    
    print(f"For the given integer k = {k}, we compute the limit.")
    print(f"The formula for the limit is (k - 1) / (k + 1).")
    print("\n--- Calculation ---")
    print(f"Numerator = k - 1 = {k} - 1 = {numerator}")
    print(f"Denominator = k + 1 = {k} + 1 = {denominator}")
    
    if denominator == 0:
        result_str = "undefined (division by zero)"
    else:
        result = numerator / denominator
        result_str = f"{result:.4f}"

    print(f"Final Equation: ({k} - 1) / ({k} + 1) = {numerator} / {denominator}")
    print(f"The value of the limit is: {result_str}")

# Example usage with a value for k.
# You can change this value to any integer >= 2.
k_value = 4
compute_limit(k_value)