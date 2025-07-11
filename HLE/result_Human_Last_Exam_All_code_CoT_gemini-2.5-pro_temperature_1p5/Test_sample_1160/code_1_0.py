import math

def calculate_limit_expression(n, C=1.0):
    """
    Calculates the expression for a given n, based on the model p_n = C/n.

    The problem asks to find the limit of log(1/p_n) / log(n) as n -> infinity.
    Our physical argument suggests that p_n, the escape probability, behaves like C/n for large n.
    This function demonstrates that for such a p_n, the limit is 1.
    """
    if n <= 1:
        print(f"n must be greater than 1.")
        return

    # Step 1: Calculate p_n based on the derived model p_n = C/n
    p_n = C / n

    # Step 2: Calculate 1/p_n
    inv_p_n = 1.0 / p_n

    # Step 3: Calculate the numerator log(1/p_n)
    log_inv_p_n = math.log(inv_p_n)

    # Step 4: Calculate the denominator log(n)
    log_n = math.log(n)

    # Step 5: Calculate the final ratio
    ratio = log_inv_p_n / log_n

    print(f"For n = {n}:")
    print(f"  p_n = C/n = {C}/{n} = {p_n:.6f}")
    print(f"  1/p_n = {inv_p_n:.2f}")
    print(f"  ln(1/p_n) = ln({inv_p_n:.2f}) = {log_inv_p_n:.6f}")
    print(f"  ln(n) = ln({n}) = {log_n:.6f}")
    print(f"  Final Ratio = {log_inv_p_n:.6f} / {log_n:.6f} = {ratio:.6f}\n")

def main():
    """
    Main function to run the demonstration.
    """
    print("This script demonstrates the calculation of the limit lim_{n->inf} [ln(1/p_n) / ln(n)].")
    print("Based on our analysis, p_n is approximately C/n, which leads to a limit of 1.")
    print("We show the calculation for several large values of n, assuming C=1 for simplicity.")
    print("-" * 70)

    # Calculate for several large n values to show the convergence to 1
    test_values = [10, 100, 1000, 10000, 100000]
    for n_val in test_values:
        calculate_limit_expression(n_val)

    print("As n gets larger, the ratio approaches 1.")

if __name__ == "__main__":
    main()
