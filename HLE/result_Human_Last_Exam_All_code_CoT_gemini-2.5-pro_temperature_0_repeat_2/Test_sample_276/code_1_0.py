import random
import math

def solve():
    """
    Calculates the probability that a circle of radius 6 thrown onto a unit grid
    intersects exactly 47 cells, using a Monte Carlo simulation.
    """
    # The number of intersected cells is N_cells = 48 - N_vertex.
    # For N_cells = 47, we need N_vertex = 1.
    # This simulation finds the probability that the circle contains exactly one vertex.

    R = 6.0
    R_squared = R**2
    num_trials = 2 * 10**6  # A large number of trials for better accuracy
    count_n_vertex_is_1 = 0

    # Define a safe range of integer coordinates (k, l) to check.
    # For a center (u,v) in [0,1]x[0,1] and R=6, the circle can only contain
    # vertices (k,l) where k and l are roughly between -6 and 6.
    k_min, k_max = -7, 8
    l_min, l_max = -7, 8

    for _ in range(num_trials):
        # Generate a random center for the circle in the unit square
        u = random.random()
        v = random.random()

        n_vertex = 0
        # Count the number of grid vertices inside the circle
        for k in range(k_min, k_max):
            for l in range(l_min, l_max):
                dist_sq = (k - u)**2 + (l - v)**2
                if dist_sq < R_squared:
                    n_vertex += 1
        
        # Check if the condition (N_vertex = 1) is met
        if n_vertex == 1:
            count_n_vertex_is_1 += 1

    probability = count_n_vertex_is_1 / num_trials

    print("The equation for the number of intersected cells (N_cells) is:")
    print("N_cells = 2 * (Number of vertical lines crossed + Number of horizontal lines crossed) - Number of grid vertices inside the circle")
    print("For a circle of radius 6, this becomes (almost everywhere):")
    print("N_cells = 2 * (12 + 12) - N_vertex = 48 - N_vertex")
    print("We want to find the probability that N_cells = 47.")
    print("47 = 48 - N_vertex")
    print("This implies the number of vertices inside the circle must be 1 (N_vertex = 1).")
    print(f"\nThe simulation with {num_trials} trials estimates the probability of this event.")
    print(f"Approximate probability: {probability:.4f}")

solve()