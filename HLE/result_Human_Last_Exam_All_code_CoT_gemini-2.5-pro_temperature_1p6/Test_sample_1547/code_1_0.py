import random

def is_in_transient(pos):
    """Check if position is in the transient set (the origin)."""
    return pos == (0, 0, 0)

def is_in_recurrent(pos):
    """Check if position is in the recurrent set (the plane z=0)."""
    return pos[2] == 0

def simulate_walk(start_pos, is_in_target, escape_norm):
    """
    Simulates a single 3D random walk.
    Returns True if it hits the target set, False if it escapes.
    An L-infinity norm (max coordinate) is used for the escape condition.
    """
    pos = list(start_pos)
    moves = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    
    # A generous max_steps limit to prevent near-infinite loops.
    # The time to travel a distance R is typically of the order R^2.
    max_steps = 4 * escape_norm * escape_norm

    for _ in range(max_steps):
        if is_in_target(tuple(pos)):
            return True
        
        if max(abs(p) for p in pos) >= escape_norm:
            return False

        move = random.choice(moves)
        pos[0] += move[0]
        pos[1] += move[1]
        pos[2] += move[2]
            
    return False # Assume escape if max_steps is reached

def estimate_hitting_prob(start_pos, is_in_target, n_walks):
    """
    Estimates the hitting probability by running many simulations.
    The escape boundary is set far from the starting point.
    """
    hits = 0
    start_norm = max(abs(p) for p in start_pos) if start_pos != (0,0,0) else 1
    escape_norm = max(20, 5 * start_norm)

    for _ in range(n_walks):
        if simulate_walk(start_pos, is_in_target, escape_norm):
            hits += 1
    
    return hits / n_walks

def main():
    """
    Main function to run the simulations and print the results.
    """
    n_walks = 5000
    distances = [1, 5, 10, 15, 20]

    print("--- Numerical Illustration of Hitting Probabilities in 3D ---")
    print(f"Estimations based on {n_walks} simulated walks per point.")
    print("-" * 60)
    print("Set A_transient = The origin {(0,0,0)}")
    print("Set A_recurrent = The plane {(x,y,0) | x,y in Z}")
    print("-" * 60)

    for dist in distances:
        start_point = (dist, 0, 0)
        
        # Calculate for the transient set
        prob_transient = estimate_hitting_prob(start_point, is_in_transient, n_walks)
        print(f"For distance = {dist:<2}:")
        print(f"  P(hit A_transient from {start_point}) ≈ {prob_transient:.4f}")

        # Calculate for the recurrent set
        # For a plane, the hitting probability is theoretically 1.
        # Simulation may slightly underestimate it due to the finite escape boundary.
        prob_recurrent = estimate_hitting_prob(start_point, is_in_recurrent, n_walks)
        print(f"  P(hit A_recurrent from {start_point}) ≈ {prob_recurrent:.4f}")
        print("-" * 60)

if __name__ == "__main__":
    main()