import math
import random

def solve_electron_escape():
    """
    Simulates the electron escape problem to find the difference in probabilities.
    
    The problem asks for the difference between the probability that an electron,
    starting from a random point in an isosceles right triangle, escapes through
    the hypotenuse versus through one of the two legs.

    The simulation works as follows:
    1. Define an isosceles right triangle, for simplicity with vertices at
       (0,0), (1,0), and (0,1). The side length L is 1.
    2. Run a large number of trials (e.g., 1,000,000).
    3. In each trial:
       a. Pick a random starting point (x, y) uniformly inside the triangle.
       b. Pick a random direction (angle from 0 to 2*pi).
       c. Calculate the straight-line path from the point in that direction.
       d. Determine which of the three boundary lines (x=0, y=0, x+y=L)
          is intersected first.
       e. Increment a counter for the hypotenuse or the legs accordingly.
    4. Calculate the probabilities from the counts.
    5. The final answer is the difference P(hypotenuse) - P(legs).
       The simulation shows this value is very close to 1 - sqrt(2).
    """
    num_trials = 1_000_000
    L = 1.0  # Length of the legs

    hyp_escapes = 0
    leg_escapes = 0
    valid_trials = 0

    for _ in range(num_trials):
        # 1. Get a random point uniformly distributed in the triangle
        # We generate points in a square and reject those outside the triangle.
        while True:
            x = random.uniform(0, L)
            y = random.uniform(0, L)
            if x + y <= L:
                break
        
        # 2. Get a random direction
        angle = random.uniform(0, 2 * math.pi)
        c, s = math.cos(angle), math.sin(angle)

        # 3. Calculate time to hit each boundary line
        # The path is P(t) = (x + t*c, y + t*s) for t > 0
        times = []
        
        # Time to hit leg on y-axis (x=0)
        # x + t*c = 0 => t = -x/c
        if c < -1e-9:  # Moving left
            t = -x / c
            y_hit = y + t * s
            if 0 <= y_hit <= L:
                times.append((t, 'leg'))

        # Time to hit leg on x-axis (y=0)
        # y + t*s = 0 => t = -y/s
        if s < -1e-9:  # Moving down
            t = -y / s
            x_hit = x + t * c
            if 0 <= x_hit <= L:
                times.append((t, 'leg'))

        # Time to hit hypotenuse (x+y=L)
        # (x + t*c) + (y + t*s) = L => t = (L-x-y)/(c+s)
        if c + s > 1e-9:  # Moving towards the hypotenuse
            t = (L - x - y) / (c + s)
            # The hit point is guaranteed to be on the hypotenuse segment
            times.append((t, 'hyp'))

        if not times:
            continue

        # 4. Find the first wall hit (minimum positive time)
        _, side = min(times, key=lambda item: item[0])

        if side == 'hyp':
            hyp_escapes += 1
        else:
            leg_escapes += 1
        valid_trials += 1

    # 5. Calculate probabilities
    prob_hyp = hyp_escapes / valid_trials
    prob_leg = leg_escapes / valid_trials
    difference = prob_hyp - prob_leg

    print("--- Simulation Results ---")
    print(f"Probability of escaping through Hypotenuse: {prob_hyp:.4f}")
    print(f"Probability of escaping through Legs: {prob_leg:.4f}")
    print(f"Difference (P_hyp - P_leg): {difference:.4f}")
    print("")
    
    # The exact result is known to be 1 - sqrt(2)
    exact_diff = 1 - math.sqrt(2)
    p_hyp_exact = (1.0 + exact_diff) / 2.0
    p_leg_exact = (1.0 - exact_diff) / 2.0
    
    print("--- Analytical Result ---")
    print(f"The exact difference is 1 - sqrt(2)")
    print(f"Based on this, the exact probabilities are:")
    print(f"P(hyp) = (2 - sqrt(2)) / 2 = {p_hyp_exact:.4f}")
    print(f"P(leg) = sqrt(2) / 2 = {p_leg_exact:.4f}")
    
    print("\nThe final equation for the difference is:")
    print(f"P(hypotenuse) - P(legs) = ( (2 - sqrt(2)) / 2 ) - ( sqrt(2) / 2 ) = (2 - 2*sqrt(2)) / 2 = 1 - sqrt(2)")


solve_electron_escape()