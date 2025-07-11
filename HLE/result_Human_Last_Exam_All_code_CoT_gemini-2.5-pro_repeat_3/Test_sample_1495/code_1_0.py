import random

def calculate_area_ratio(d, e, f):
    """
    Calculates the ratio Area(XYZ)/Area(ABC) based on a simplified
    formula derived from Routh's theorem.

    Args:
        d: A float in [0,1], representing the ratio BD/BC.
        e: A float in [0,1], representing the ratio CE/CA.
        f: A float in [0,1], representing the ratio AF/AB.

    Returns:
        The ratio of the areas, a float.
    """
    # This formula is a stable version of Routh's Theorem's formula.
    # Numerator of the ratio:
    term1 = d * e * f
    term2 = (1 - d) * (1 - e) * (1 - f)
    numerator = (term1 - term2)**2

    # Denominator of the ratio:
    den1 = 1 - d + d * e
    den2 = 1 - e + e * f
    den3 = 1 - f + f * d

    # The denominator can only be zero in edge cases with probability 0
    # for continuous random variables (e.g., d=1 and e=0).
    # We handle this to prevent division by zero, though it's unlikely to occur.
    denominator = den1 * den2 * den3
    if denominator == 0:
        return 0.0

    return numerator / denominator

def monte_carlo_simulation(num_samples):
    """
    Runs a Monte Carlo simulation to estimate the desired probability.

    Args:
        num_samples: The number of random trials to perform.

    Returns:
        The estimated probability.
    """
    total_ratio_sum = 0
    for _ in range(num_samples):
        # Generate three independent random numbers uniformly in [0, 1]
        d = random.random()
        e = random.random()
        f = random.random()
        
        # Calculate the area ratio for this trial and add to the sum
        total_ratio_sum += calculate_area_ratio(d, e, f)

    # The estimated probability is the average of all calculated ratios
    return total_ratio_sum / num_samples

if __name__ == "__main__":
    # Number of samples for the simulation. More samples lead to a more
    # accurate estimate but take longer to run.
    num_samples = 2000000

    # Run the simulation
    estimated_probability = monte_carlo_simulation(num_samples)

    print(f"Based on {num_samples} samples, the estimated probability is: {estimated_probability}")
    
    # The problem has a known analytical solution, which we can state for comparison.
    theoretical_probability = 1/10
    print(f"The exact theoretical probability is 1/10 = {theoretical_probability}")
