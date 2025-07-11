import random
import math

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D random walk visits a target set
    before leaving a specified disk, using a Monte Carlo simulation.
    """
    # --- Simulation Parameters ---
    # Number of random walks to simulate for statistical accuracy.
    num_simulations = 2_000_000

    # Starting position of the random walk.
    start_pos = (0, 300)

    # Target set of points to visit. Using a set for fast lookups.
    target_set = {(0, 0), (2, 0)}

    # The walk stops if it leaves the disk of radius 1000.
    # We check against the radius squared to avoid costly square root operations.
    radius_squared = 1000**2

    # --- Simulation Logic ---
    successful_hits = 0

    # Possible moves (dx, dy) for the walker.
    moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    # Main simulation loop.
    for _ in range(num_simulations):
        # Start each walk at the specified starting position.
        x, y = start_pos

        # Run the walk until it hits the target or leaves the disk.
        while True:
            # Check for success: Has the walker reached the target?
            if (x, y) in target_set:
                successful_hits += 1
                break  # End this walk as a success.

            # Check for failure: Has the walker left the disk?
            if x*x + y*y > radius_squared:
                break  # End this walk as a failure.

            # If not terminated, take one random step.
            dx, dy = random.choice(moves)
            x += dx
            y += dy

    # --- Calculate and Print Result ---
    # The probability is the ratio of successful walks to the total number of walks.
    probability = successful_hits / num_simulations

    print(f"Simulation Details:")
    print(f"Starting Position: {start_pos}")
    print(f"Target Set: {target_set}")
    print(f"Boundary: Disk of radius {int(math.sqrt(radius_squared))}")
    print("-" * 20)
    
    # Outputting the numbers from the final calculation as requested.
    print(f"Final Calculation:")
    print(f"Successful walks: {successful_hits}")
    print(f"Total walks: {num_simulations}")
    print(f"Probability = {successful_hits} / {num_simulations} = {probability}")
    print("-" * 20)

    # Format the final answer to three significant digits.
    if probability > 0:
        sig_digits = 3
        # This formula calculates the number of decimal places needed for a given
        # number of significant figures.
        num_decimal_places = -int(math.floor(math.log10(abs(probability)))) + (sig_digits - 1)
        formatted_prob = f"{probability:.{num_decimal_places}f}"
    else:
        formatted_prob = "0.00"

    print(f"The probability rounded to three significant digits is: {formatted_prob}")


if __name__ == "__main__":
    solve_random_walk_probability()
