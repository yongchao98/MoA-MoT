import math
import random

def solve():
    """
    Calculates the probability that a circumference of radius R=6 thrown randomly
    onto a grid of unit squares intersects exactly 47 squares, using a Monte Carlo simulation.
    """
    R = 6.0
    R2 = R * R
    N_TARGET = 47
    
    # We can model the random throw by placing the circle's center (x, y)
    # uniformly in the unit square [0, 1] x [0, 1].
    
    # A large number of trials is needed for precision. 10 million trials
    # should provide sufficient stability for four significant digits.
    n_trials = 10_000_000
    count_target = 0

    # The circle with center (x,y) in [0,1]x[0,1] and radius 6
    # can only intersect squares [k, k+1]x[l, l+1] where k and l are in a limited range.
    # The circle's x-extent is [x-R, x+R] which is within [-6, 7].
    # So a square can only be intersected if its k-index is between -7 and 7. Same for l.
    k_min, k_max = -7, 7
    l_min, l_max = -7, 7

    for _ in range(n_trials):
        # 1. Generate a random center (x, y)
        x = random.random()
        y = random.random()

        # 2. Count the number of intersected squares
        intersected_squares_count = 0
        
        for k in range(k_min, k_max + 1):
            for l in range(l_min, l_max + 1):
                # A square [k, k+1] x [l, l+1] is intersected by the circumference if
                # min_dist(center, square) <= R <= max_dist(center, square).

                # Find the point (px, py) in the square closest to the center (x, y)
                px = max(k, min(k + 1, x))
                py = max(l, min(l + 1, y))
                
                min_dist_sq = (px - x)**2 + (py - y)**2
                
                # If the circle is completely outside the square, no intersection.
                if min_dist_sq > R2:
                    continue
                
                # Find the vertex of the square furthest from the center (x, y)
                # to determine the maximum distance.
                furthest_x = k if x > k + 0.5 else k + 1
                furthest_y = l if y > l + 0.5 else l + 1

                max_dist_sq = (furthest_x - x)**2 + (furthest_y - y)**2
                
                # If the circle completely contains the square, the circumference doesn't intersect it.
                if max_dist_sq < R2:
                    continue

                # Otherwise, the circumference must cross the square.
                intersected_squares_count += 1
        
        if intersected_squares_count == N_TARGET:
            count_target += 1
            
    probability = count_target / n_trials
    
    print("Plan execution:")
    print(f"Radius of circumference (R): {R}")
    print(f"Target number of intersected cells: {N_TARGET}")
    print(f"Total number of random throws (trials): {n_trials}")
    print(f"Number of throws where the circle intersected {N_TARGET} cells: {count_target}")
    print(f"Probability = {count_target} / {n_trials}")

    # Format the final result to four significant digits.
    final_answer_str = f"{probability:.4g}"
    print(f"\nThe approximate probability is: {final_answer_str}")
    
    # Do not remove the following line
    print(f"\n<<<{final_answer_str}>>>")

solve()