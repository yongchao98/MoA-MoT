import math

def solve_for_p():
    """
    Solves the equation p^3 + p - 1 = 0 using Newton's method.
    f(p) = p^3 + p - 1
    f'(p) = 3p^2 + 1
    Iteration: p_new = p_old - f(p_old) / f'(p_old)
    """
    # Start with an initial guess
    p = 0.5
    # Iterate 10 times for high precision
    for _ in range(10):
        p_cubed = p * p * p
        p = p - (p_cubed + p - 1) / (3 * p * p + 1)
    return p

def main():
    """
    Calculates and prints the minimal cost coefficient.
    """
    # Step 1: Find the optimal probability 'p' by solving p^3 + p - 1 = 0
    p_optimal = solve_for_p()

    # Step 2: Calculate the minimal cost per bit of information
    # The formula is C = -1 / log2(p)
    # This is equivalent to C = -ln(2) / ln(p)
    log_base2_of_p = math.log2(p_optimal)
    min_cost_per_bit = -1 / log_base2_of_p

    # Step 3: Print the results as requested.
    # The "final equation" is the one for the cost. We will show the numbers involved.
    print(f"The optimal 'yes' probability 'p' is the root of the equation p^3 + p - 1 = 0.")
    print(f"Solving numerically, we find p = {p_optimal:.5f}")
    print("\nThe minimal cost per bit of information is given by the equation: C = -1 / log2(p)")
    print("Substituting the value of p:")
    print(f"C = -1 / log2({p_optimal:.5f})")
    print(f"C = -1 / {log_base2_of_p:.5f}")
    print(f"C = {min_cost_per_bit:.5f}")
    print("\nRounding to 3 decimal places, the minimal cost coefficient is:")
    print(f"{min_cost_per_bit:.3f}")

if __name__ == "__main__":
    main()