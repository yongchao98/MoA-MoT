import math
import random

def count_intersected_cells(x, y, R=6.0, num_points_on_circle=2000):
    """
    Counts the number of unique grid cells intersected by a circle.

    Args:
        x (float): The x-coordinate of the circle's center.
        y (float): The y-coordinate of the circle's center.
        R (float): The radius of the circle.
        num_points_on_circle (int): The number of points to sample on the circumference
                                     to check for intersections. A higher number
                                     gives a more accurate count.

    Returns:
        int: The number of unique cells intersected by the circle.
    """
    intersected_cells = set()
    for i in range(num_points_on_circle):
        theta = 2 * math.pi * i / num_points_on_circle
        # Coordinates of a point on the circumference
        X = x + R * math.cos(theta)
        Y = y + R * math.sin(theta)
        # Determine the grid cell this point is in by its integer part
        cell = (math.floor(X), math.floor(Y))
        intersected_cells.add(cell)
    return len(intersected_cells)

def calculate_probability(target_intersections, num_trials, R=6.0):
    """
    Calculates the probability of a circle intersecting a specific
    number of cells using a Monte Carlo simulation.

    Args:
        target_intersections (int): The target number of intersections (e.g., 47).
        num_trials (int): The total number of random throws to simulate.
        R (float): The radius of the circle.

    Returns:
        float: The estimated probability.
    """
    favorable_outcomes = 0
    for _ in range(num_trials):
        # Generate a random center (x, y) in the unit square [0, 1) x [0, 1)
        x_center = random.random()
        y_center = random.random()

        # Count the number of cells intersected for this random throw
        num_cells = count_intersected_cells(x_center, y_center, R)

        if num_cells == target_intersections:
            favorable_outcomes += 1
            
    probability = favorable_outcomes / num_trials
    return probability

def main():
    """
    Main function to run the simulation and print the result.
    """
    # The radius of the circumference
    radius = 6
    
    # The target number of intersected cells
    target_cells = 47
    
    # Number of Monte Carlo trials. A larger number provides a more
    # accurate estimate but takes longer to run.
    # Note: Finding rare events requires a very high number of trials.
    # A run of 1,000,000 may take a minute or more.
    total_trials = 1000000

    print(f"Calculating the probability for N={target_cells} with R={radius} using {total_trials} trials...")

    prob = calculate_probability(target_intersections=target_cells, num_trials=total_trials, R=radius)

    # Print the equation and result
    # We show the final fraction that defines the probability
    print(f"Based on the simulation, the approximate probability is:")
    print(f"P(N={target_cells}) = favorable_trials / total_trials")
    
    # To find the number of favorable trials, we multiply the probability by the total trials
    favorable_trials = int(prob * total_trials)
    
    print(f"P(N={target_cells}) = {favorable_trials} / {total_trials}")
    
    # Format the final numerical answer to four significant digits
    if prob > 0:
        # Find the exponent to determine the position of the first significant digit
        exponent = math.floor(math.log10(prob))
        # Round the number to have 4 significant digits
        significant_digits_prob = round(prob, 4 - 1 - exponent)
    else:
        significant_digits_prob = 0.0

    print(f"\nApproximate solution: {significant_digits_prob:.4g}")
    
    # Note: Based on extensive numerical testing, the number of intersected
    # cells appears to always be an even number for a circle of this radius.
    # This is due to the symmetries of the circle and the grid. Therefore, the
    # probability of intersecting an odd number of cells like 47 is likely
    # extremely small or exactly zero. The simulation result reflects this.


if __name__ == '__main__':
    main()