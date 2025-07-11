import math

def calculate_lowest_coverage(n, alpha):
    """
    Calculates the lowest guaranteed coverage probability for a LOO-based conformal interval.

    Args:
        n (int): The number of training data points.
        alpha (float): The significance level, e.g., 0.1 for 90% nominal coverage.

    Returns:
        float: The lowest possible coverage probability.
    """
    # The number of scores in total (n training scores + 1 test score)
    num_total_points = n + 1

    # The formula for the lowest guaranteed coverage is 1 - floor((n+1)*alpha) / (n+1)
    numerator = math.floor(num_total_points * alpha)
    denominator = num_total_points
    
    miscoverage_prob = numerator / denominator
    coverage_prob = 1 - miscoverage_prob

    print(f"For n = {n} and alpha = {alpha}:")
    print(f"The total number of points considered is n + 1 = {num_total_points}")
    print(f"The number of allowed miscoverages is floor((n+1) * alpha) = floor({num_total_points} * {alpha}) = floor({num_total_points * alpha}) = {numerator}")
    print(f"The maximum miscoverage probability is {numerator}/{denominator} = {miscoverage_prob}")
    print("The lowest coverage probability is given by the formula: 1 - floor((n+1)*alpha) / (n+1)")
    print(f"Final Calculation: 1 - {numerator} / {denominator} = {coverage_prob}")
    print("\n---")
    print("The final expression for the lowest possible coverage probability is:")
    print("1 - floor((n + 1) * alpha) / (n + 1)")


# Example from the text: n=19, alpha=0.1
calculate_lowest_coverage(19, 0.1)

# Example from the text: n=20, alpha=0.1
calculate_lowest_coverage(20, 0.1)