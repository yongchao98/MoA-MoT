import math

def calculate_lowest_coverage(n: int, alpha: float):
    """
    Calculates the guaranteed minimum coverage probability for a CV+ prediction interval.

    The formula for the sharp lower bound on the coverage probability is:
    P(coverage) >= floor((n + 1) * (1 - alpha)) / (n + 1)

    Args:
        n (int): The number of training samples.
        alpha (float): The desired significance level (e.g., 0.1 for 90% intervals).
    """
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return
    if not (0 < alpha < 1):
        print("Error: alpha must be between 0 and 1.")
        return

    # Numerator of the formula
    numerator_val = math.floor((n + 1) * (1 - alpha))
    
    # Denominator of the formula
    denominator_val = n + 1
    
    # Calculate the lowest possible coverage probability
    lowest_coverage = numerator_val / denominator_val

    print("The method described is known as CV+ (or Jackknife+).")
    print("For any data distribution and any algorithm, there is a guaranteed lower bound on the coverage probability.")
    print("\nThe formula for this lowest possible probability is:")
    print("floor((n + 1) * (1 - alpha)) / (n + 1)\n")

    print(f"Given the values:")
    print(f"n = {n}")
    print(f"alpha = {alpha}\n")

    print("The calculation is:")
    # Using 'math.floor' in the string to represent the floor operation
    print(f"floor(({n} + 1) * (1 - {alpha})) / ({n} + 1)")
    print(f"= floor({n + 1} * {1 - alpha}) / {n + 1}")
    print(f"= floor({(n + 1) * (1 - alpha)}) / {n + 1}")
    print(f"= {numerator_val} / {denominator_val}")
    print(f"= {lowest_coverage:.4f}\n")
    
    print(f"Therefore, the lowest possible value for the coverage probability is {lowest_coverage:.4f}, or {lowest_coverage:%}.")

# --- Example Usage ---
# You can change these values to see the result for different scenarios.
number_of_samples = 100
significance_level = 0.1

calculate_lowest_coverage(number_of_samples, significance_level)