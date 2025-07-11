import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the derived formula.

    The problem involves a maximum of N random variables, where N itself is a random variable.
    - Magnitudes X ~ Pareto(2), so CDF F_X(x) = 1 - 1/x^2.
    - Number of years N ~ LogSeries(1/2).
    - The expected maximum M = max(X_1, ..., X_N) is found to be E[M] = pi / (2 * log(2)) - 1.
    This script calculates and prints the value based on this final formula.
    """
    # --- 1. Define constants ---
    pi_val = math.pi
    log2_val = math.log(2)

    # --- 2. Calculate the components of the formula ---
    denominator = 2 * log2_val
    quotient = pi_val / denominator
    result = quotient - 1

    # --- 3. Print the derivation step-by-step as requested ---
    print("The analytical solution for the expected maximum magnitude, E[M], is:")
    print("E[M] = pi / (2 * log(2)) - 1")
    print("\nSubstituting the values of the constants:")
    print(f"E[M] = {pi_val} / (2 * {log2_val}) - 1")
    print(f"E[M] = {pi_val} / {denominator} - 1")
    print(f"E[M] = {quotient} - 1")
    print(f"\nThe final expected maximum magnitude is:")
    print(f"E[M] = {result}")

if __name__ == "__main__":
    solve_earthquake_magnitude()