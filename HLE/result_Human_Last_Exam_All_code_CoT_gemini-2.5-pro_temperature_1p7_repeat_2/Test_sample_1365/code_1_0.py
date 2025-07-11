import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes for the given
    variant of the experts problem.

    The formula for the bound is M <= (c - 1) + 2 * c * ln((n + 1) / 2).

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes after which an expert is removed.

    Returns:
        float: The calculated upper bound for the number of mistakes.
    """
    if n <= 0 or c <= 1:
        # The model assumes n >= 1 and c > 1. The true expert makes < c mistakes.
        # If c=1, true expert makes 0 mistakes. Any other expert is removed on 1st mistake.
        # If n=1, the algorithm is just the true expert, making < c mistakes.
        if n == 1:
            return float(c - 1)
        else:
            raise ValueError("n must be positive, and c must be greater than 1.")

    # M1 is bounded by the number of mistakes of the true expert
    m1_bound = c - 1
    
    # M2 is bounded by 2 * c * ln((n + 1) / 2)
    # The argument to log must be > 0. n>=1 ensures (n+1)/2 >= 1.
    log_term = math.log((n + 1) / 2)
    m2_bound = 2 * c * log_term

    total_bound = m1_bound + m2_bound
    
    return total_bound

def main():
    """
    Main function to demonstrate the calculation of the mistake bound.
    """
    # Example values for n and c
    n = 21
    c = 10

    print(f"Calculating the mistake bound for n = {n} experts and c = {c} mistakes for removal.")
    print("-" * 30)
    
    # Breakdown of the formula
    print("The upper bound M is given by the formula:")
    print("M <= M1 + M2")
    print("where M1 is the bound for mistakes when the true expert is wrong,")
    print("and M2 is the bound for mistakes when the true expert is right.")
    print("\nStep 1: Bound M1")
    print(f"M1 <= c - 1")
    print(f"M1 <= {c} - 1 = {c - 1}")

    print("\nStep 2: Bound M2")
    print("M2 <= 2 * c * ln((n + 1) / 2)")
    log_val = math.log((n + 1) / 2)
    print(f"M2 <= 2 * {c} * ln(({n} + 1) / 2)")
    print(f"M2 <= {2 * c} * ln({(n + 1) / 2})")
    print(f"M2 <= {2 * c} * {log_val:.4f} = {(2 * c * log_val):.4f}")
    
    # Calculate final bound
    bound = calculate_mistake_bound(n, c)

    print("\nStep 3: Combine bounds for the final result")
    print("M <= (c - 1) + 2 * c * ln((n + 1) / 2)")
    print(f"M <= ({c} - 1) + 2 * {c} * ln(({n} + 1) / 2)")
    print(f"M <= {c - 1} + {(2 * c * log_val):.4f}")
    print(f"M <= {bound:.4f}")

    print("\nFinal calculated upper bound on mistakes:")
    print(bound)


if __name__ == "__main__":
    main()
