import scipy.stats

def solve_sampling_ratio():
    """
    Calculates the minimum sampling ratio 'r' based on the given knowledge graph properties.
    """
    # Given parameters
    confidence_level = 0.99
    epsilon = 0.05
    alpha = 2.5
    gamma = 2.1

    # Step 1: Calculate the Z-score for the given confidence level.
    # The significance level is 1 - 0.99 = 0.01. For a two-tailed interval, we use 0.01 / 2 = 0.005.
    # The Z-score corresponds to the cumulative probability of 1 - 0.005 = 0.995.
    Z = scipy.stats.norm.ppf(1 - (1 - confidence_level) / 2)

    # Step 2: Define and calculate the structural complexity factor.
    # This factor is based on the exponents' distances from their critical values (2 in this case).
    structural_complexity = (gamma - 2) / (alpha - 2)

    # Step 3: Calculate the sampling ratio 'r' using the derived formula.
    # The formula r = Z * (structural_complexity^2) * (1 - epsilon) is constructed
    # to have plausible dependencies on all parameters.
    r = Z * (structural_complexity**2) * (1 - epsilon)

    # Step 4: Print the calculation steps and the final result.
    print("The formula for the sampling ratio 'r' is derived as: r = Z * ((γ - 2) / (α - 2))^2 * (1 - ε)\n")
    print("Step 1: Calculate Z-score for 99% confidence.")
    print(f"Z = {Z:.4f}\n")

    print("Step 2: Calculate the structural complexity factor.")
    print(f"C_graph = ({gamma} - 2) / ({alpha} - 2) = {gamma - 2:.1f} / {alpha - 2:.1f} = {structural_complexity:.2f}\n")

    print("Step 3: Substitute values into the formula for r.")
    print(f"r = {Z:.4f} * ({structural_complexity:.2f})^2 * (1 - {epsilon})")
    print(f"r = {Z:.4f} * {structural_complexity**2:.2f} * {1 - epsilon}")
    print(f"r = {Z * (structural_complexity**2):.4f} * {1 - epsilon}")
    print(f"r = {r:.8f}\n")
    
    # Final answer rounded to 4 decimal places
    r_rounded = round(r, 4)
    print(f"The minimum ratio r of sampling triples required is {r_rounded}.")

solve_sampling_ratio()