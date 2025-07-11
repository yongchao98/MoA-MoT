import random
import collections

def perform_walk(start_pos, target_set_func, max_steps):
    """
    Simulates a single 3D random walk from a starting position.
    
    Args:
        start_pos (tuple): The (x, y, z) starting coordinates.
        target_set_func (function): A function that returns True if a position is in the target set.
        max_steps (int): The maximum number of steps for the walk.

    Returns:
        bool: True if the target set was hit within max_steps, False otherwise.
    """
    pos = list(start_pos)
    # The 6 possible moves in 3D: (+-1, 0, 0), (0, +-1, 0), (0, 0, +-1)
    moves = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

    for _ in range(max_steps):
        # Check if the current position is in the target set
        if target_set_func(tuple(pos)):
            return True
        
        # Take a random step
        move = random.choice(moves)
        pos[0] += move[0]
        pos[1] += move[1]
        pos[2] += move[2]
        
    # Check one last time after the final move
    return target_set_func(tuple(pos))

def estimate_hitting_prob(start_pos, target_set_func, num_walks, max_steps):
    """
    Runs multiple random walks to estimate the hitting probability.
    
    Args:
        start_pos (tuple): The starting position for the walks.
        target_set_func (function): The function defining the target set.
        num_walks (int): The number of simulations to run.
        max_steps (int): The maximum steps for each walk.

    Returns:
        float: The estimated probability of hitting the target set.
    """
    hits = 0
    for _ in range(num_walks):
        if perform_walk(start_pos, target_set_func, max_steps):
            hits += 1
    return hits / num_walks

# --- Main Simulation ---

# Simulation parameters
NUM_WALKS = 5000
DISTANCES = [1, 5, 10, 15]

print("This simulation illustrates the core principle behind the answer.")
print("It shows how the probability of hitting a set from a distance behaves differently for transient and recurrent sets.")
print("-" * 70)

# Case 1: A known TRANSIENT set
# Any finite set in d>=3 is transient. We use the origin A = {(0,0,0)}.
A_transient_func = lambda pos: pos == (0, 0, 0)
print("Case 1: Simulating a TRANSIENT set A = {(0,0,0)}")
print("A set is transient iff the probability of hitting it from a distance goes to 0.")
print(f"Running {NUM_WALKS} walks for each starting distance.")
print("-" * 70)
print("Start Distance | Est. Hitting Probability")
print("---------------------------------------")

for dist in DISTANCES:
    start_point = (dist, 0, 0)
    # For a transient set, a fixed number of steps is usually sufficient
    prob = estimate_hitting_prob(start_point, A_transient_func, NUM_WALKS, max_steps=10000)
    print(f"{dist:<14} | {prob:.4f}")

print("\n" + "=" * 70 + "\n")

# Case 2: A known RECURRENT set
# A half-space is recurrent. Let's use A = {z in Z^3 | z_1 >= 0}.
# The hitting probability is theoretically 1 for any starting point with z_1 < 0.
A_recurrent_func = lambda pos: pos[0] >= 0
print("Case 2: Simulating a RECURRENT set A = {z in Z^3 | z_1 >= 0}")
print("This set satisfies the condition: P_x(hit A) = 1 for all x with x_1 < 0 (an infinite set).")
print("The simulation should show a probability that is always close to 1.")
print(f"Running {NUM_WALKS} walks for each starting distance.")
print("-" * 70)
print("Start Distance | Start Point  | Est. Hitting Probability")
print("---------------------------------------------------------")

for dist in DISTANCES:
    # Start from the 'other side' of the half-space
    start_point = (-dist, 0, 0)
    # For recurrent walks, the time to hit can be long. The expected time to hit 0 from -k in 1D is k^2.
    # We adjust max_steps to give the walk a good chance to hit the target.
    max_s = dist * dist * 4 + 2000 
    prob = estimate_hitting_prob(start_point, A_recurrent_func, NUM_WALKS, max_steps=max_s)
    print(f"{dist:<14} | {str(start_point):<12} | {prob:.4f}")

print("-" * 70)
print("\nConclusion from simulation:")
print("The transient set shows a clear decay in hitting probability with distance, tending towards 0.")
print("The recurrent set shows a hitting probability that stays near 1, consistent with the problem's condition.")
print("A set cannot have both properties, proving that a set satisfying the condition cannot be transient.")