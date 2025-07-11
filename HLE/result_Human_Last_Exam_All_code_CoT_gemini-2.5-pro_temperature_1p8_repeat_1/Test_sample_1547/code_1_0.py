import random
import numpy as np

def simulate_hitting_probability(d, start_pos, target_set, num_simulations, max_walk_length, escape_radius):
    """
    Estimates the probability of a d-dimensional SRW hitting a target set.

    Args:
        d (int): The dimension of the space.
        start_pos (tuple): The starting position of the walk.
        target_set (set): A set of tuples representing the target points.
        num_simulations (int): The number of walk simulations to run.
        max_walk_length (int): The maximum number of steps for each walk.
        escape_radius (float): If the walk exceeds this distance, it's considered escaped.

    Returns:
        float: The estimated hitting probability.
    """
    hits = 0
    # Generate the set of possible steps (unit vectors)
    steps = []
    for i in range(d):
        steps.append(tuple(1 if j == i else 0 for j in range(d)))
        steps.append(tuple(-1 if j == i else 0 for j in range(d)))

    for _ in range(num_simulations):
        current_pos = np.array(start_pos)
        
        # A quick check if start position is already in the target
        # For tau_A = min{n >= 1: ...}, a walk must take at least one step.
        # But if we define it as min{n >= 0: ...}, this check would be important.
        # Let's assume the question implies start_pos is not in A.
        
        for _ in range(max_walk_length):
            # Take a random step
            step = np.array(random.choice(steps))
            current_pos += step

            # Check if the new position is in the target set
            if tuple(current_pos) in target_set:
                hits += 1
                break
            
            # Check if the walk has escaped
            if np.linalg.norm(current_pos) > escape_radius:
                break
    
    return hits / num_simulations

def main():
    """
    Main function to run the simulation and print results.
    """
    d = 3  # Dimension (must be >= 3 for the walk to be transient)
    
    # A simple transient set: the origin. Any finite set would also be transient.
    target_set = {(0, 0, 0)}
    
    # Simulation parameters
    num_simulations = 20000 # Increase for more accuracy
    
    print(f"Simulating hitting probability for the transient set A = {target_set} in Z^{d}")
    print("-------------------------------------------------------------------------")
    print("The property that P_x(tau_A < inf) -> 0 as |x| -> inf is key.")
    print("This simulation demonstrates that property, making it impossible for")
    print("the hitting probability to be 1 for infinitely many points x.\n")

    # We will test for starting points with increasing distance
    for dist in range(1, 11):
        # Start at (dist, 0, 0)
        start_pos = tuple([0] * d)
        start_pos = (dist,) + start_pos[1:]
        
        # Heuristics for simulation parameters
        # To reach the origin from distance 'dist', the walk needs at least dist steps.
        # A random walk of length N travels about sqrt(N) distance.
        # So we set max_walk_length proportional to dist^2.
        max_walk_length = 50 * dist * dist
        # Set an escape radius far from both the start and the target.
        escape_radius = 2 * dist + 10

        prob = simulate_hitting_probability(
            d,
            start_pos,
            target_set,
            num_simulations,
            max_walk_length,
            escape_radius
        )
        
        # The analytical value for hitting the origin from a neighbor (dist=1)
        # is 1 - 1/G(0) where G(0) for Z^3 is ~1.516. So P ~ 0.34.
        # For large dist, P ~ 1 / ( (2d-2) * G(0) * dist ). It should decay like 1/dist.
        
        print(f"Starting Point x = {start_pos}, Distance |x| = {dist}")
        print(f"Estimated P_x(tau_A < infinity) = {prob:.4f}")
        print("-" * 20)

if __name__ == "__main__":
    main()
