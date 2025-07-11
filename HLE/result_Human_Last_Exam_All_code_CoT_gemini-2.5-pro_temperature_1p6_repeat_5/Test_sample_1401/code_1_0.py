import math

def calculate_query_lower_bound(d, c=0.5):
    """
    Calculates and explains the quasi-polynomial lower bound on the number of
    queries for learning a two-hidden-layer ReLU network with an SQ algorithm.

    The theoretical lower bound is d^{Ω(log d)}, which we represent as d^(c * log(d))
    for some constant c > 0.

    Args:
        d (int): The input dimension of the learning problem.
        c (float): A constant > 0 that depends on the specific polynomial
                   factors (network size, error, tolerance) of the problem.
                   A default value is provided for demonstration.
    """
    if not isinstance(d, int) or d <= 1:
        print("Error: The dimension 'd' must be an integer greater than 1.")
        return
    if not isinstance(c, (int, float)) or c <= 0:
        print("Error: The constant 'c' must be a positive number.")
        return

    # Calculate the components of the formula d^(c * log(d))
    log_d = math.log(d)  # math.log is the natural logarithm (ln)
    exponent = c * log_d
    # The result can be extremely large, so we handle potential overflows
    try:
        query_lower_bound = d ** exponent
        result_str = f"{query_lower_bound:,.2f}"
    except OverflowError:
        result_str = "a number too large to represent"

    # Print the final equation with each number clearly shown
    print("Based on theoretical results in statistical query (SQ) learning:")
    print("The minimum number of queries needed is quasi-polynomial in the dimension d.")
    print("-" * 60)
    print(f"Formula for the lower bound: Q = d^(c * log(d))")
    print(f"For dimension d = {d} and a constant c = {c}:")
    print("\n--- Equation Breakdown ---")
    # Output each number/component in the final equation
    print(f"1. Base (d): {d}")
    print(f"2. Constant (c): {c}")
    print(f"3. Logarithm (log(d)): {log_d:.4f}")
    print(f"4. Exponent (c * log(d)): {exponent:.4f}")
    print("\n--- Final Result ---")
    print(f"Minimum Queries (Q) >= {d}^({c} * {log_d:.4f}) = {result_str}")
    print("-" * 60)


if __name__ == '__main__':
    # You can change the dimension 'd' to see how the lower bound scales.
    example_dimension = 50
    calculate_query_lower_bound(d=example_dimension)
