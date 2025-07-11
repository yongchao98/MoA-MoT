import random
import math

def solve_random_walk_problem():
    """
    This script calculates the probability that a 2D simple random walk visits
    a target set before leaving a disk via Monte Carlo simulation.
    """

    # --- Simulation Parameters ---
    start_pos = (0, 300)
    # The target set A = {(0,0), (2,0)}
    target_set = {(0, 0), (2, 0)}
    # The walk stops if it leaves the disk of radius 1000.
    # We check against the radius squared to avoid costly square root operations.
    radius_sq = 1000**2
    # Number of simulations to run for an accurate estimate.
    # 1,000,000 runs should provide a result accurate to ~3 significant digits.
    num_simulations = 1000000

    # --- Simulation Execution ---
    hit_count = 0
    # Pre-defined moves for a simple random walk
    moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    for i in range(num_simulations):
        # Initialize the walk's position for each new simulation
        x, y = start_pos

        # A single random walk
        while True:
            # Condition 1: Check if the walk has hit the target set
            if (x, y) in target_set:
                hit_count += 1
                break  # End this walk as a "hit"

            # Condition 2: Check if the walk has left the disk
            if x*x + y*y >= radius_sq:
                break  # End this walk as an "escape"

            # If no stopping condition is met, take another random step
            dx, dy = random.choice(moves)
            x += dx
            y += dy

    # --- Calculate and Print Results ---
    probability = hit_count / num_simulations

    print("--- Monte Carlo Simulation Results ---")
    print(f"Starting point: {start_pos}")
    print(f"Target set: {target_set}")
    print(f"Escape radius: {int(math.sqrt(radius_sq))}")
    print("-" * 36)
    print("The final probability is estimated by the equation:")
    print("P = (Number of Hits) / (Total Simulations)")
    print("\nBased on the simulation:")
    print(f"Total simulations run: {num_simulations}")
    print(f"Number of times the target was hit: {hit_count}")
    print(f"The estimated probability is {hit_count} / {num_simulations} = {probability}")
    print(f"\nThe probability, to three significant digits, is {probability:.3f}")


if __name__ == "__main__":
    # Note: This script may take several minutes to run due to the large number of simulations.
    solve_random_walk_problem()
