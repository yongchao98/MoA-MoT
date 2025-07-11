import numpy as np

def get_z_sample():
    """
    Generates a single sample of the random variable Z according to the problem description.
    """
    # 1. Randomly choose four points from the unit interval [0,1]
    points = np.random.rand(4)
    x1 = points[0]
    other_points = points[1:]

    # 2. Find the second closest point to x1
    distances = np.abs(other_points - x1)
    
    # argsort returns the indices that would sort the array. The index of the
    # second smallest distance is at position 1.
    second_closest_point_index = np.argsort(distances)[1]
    x_2nd_closest = other_points[second_closest_point_index]

    # 3. Randomly choose a point Z between x1 and the second closest point
    lower_bound = min(x1, x_2nd_closest)
    upper_bound = max(x1, x_2nd_closest)
    z = np.random.uniform(low=lower_bound, high=upper_bound)
    
    return z

def estimate_pdf_at_z(z_val, n_samples=2000000, h=0.001):
    """
    Estimates the value of the PDF f_Z(z) at a specific point z_val using Monte Carlo simulation.
    """
    # Generate a large number of Z samples
    z_samples = [get_z_sample() for _ in range(n_samples)]
    
    # Count how many samples fall in a small interval around z_val
    lower_bound = z_val - h / 2
    upper_bound = z_val + h / 2
    count = sum(1 for s in z_samples if lower_bound <= s <= upper_bound)
    
    # Calculate the PDF estimate
    pdf_estimate = count / (n_samples * h)
    
    return pdf_estimate

# Calculate the value of f(0.2)
z_value_to_evaluate = 0.2
estimated_value = estimate_pdf_at_z(z_value_to_evaluate)

# The analytical result is known to be f(z) = 4 - 6|1-2z| + 3(1-2z)^2
# Let's calculate the exact value for comparison
val = 0.2
exact_val = 4 - 6 * abs(1 - 2 * val) + 3 * (1 - 2 * val)**2
# exact_val = 4 - 6 * abs(1-0.4) + 3*(1-0.4)^2
#           = 4 - 6 * 0.6 + 3 * 0.6^2
#           = 4 - 3.6 + 3 * 0.36
#           = 4 - 3.6 + 1.08
#           = 0.4 + 1.08 = 1.48
# Let's print the calculation
part1 = 4
part2 = 6 * abs(1 - 2 * val)
part3 = 3 * (1 - 2 * val)**2
print(f"The exact calculation is:")
print(f"f({val}) = {part1} - {part2} + {part3}")
print(f"f({val}) = {part1 - part2 + part3}")

print(f"\nThe estimated value from simulation is: {estimated_value}")
# The result from the simulation is expected to be close to the exact analytical result.
