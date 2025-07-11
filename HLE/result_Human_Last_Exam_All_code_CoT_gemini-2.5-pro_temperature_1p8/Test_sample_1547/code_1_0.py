import numpy as np

def simulate_hitting_probability(d, A, start_pos, num_walks, max_steps, escape_radius):
    """
    Simulates random walks to estimate the hitting probability of a set A.

    Args:
        d (int): Dimension of the space.
        A (set): The target set of points (tuples).
        start_pos (tuple): The starting position of the walk.
        num_walks (int): The number of random walks to simulate.
        max_steps (int): The maximum number of steps for each walk.
        escape_radius (float): The distance from origin to be considered "escaped".

    Returns:
        float: The estimated hitting probability.
    """
    hits = 0
    
    # Pre-compute the 2*d possible unit steps
    steps = []
    for i in range(d):
        step_plus = np.zeros(d, dtype=int)
        step_plus[i] = 1
        steps.append(step_plus)
        step_minus = np.zeros(d, dtype=int)
        step_minus[i] = -1
        steps.append(step_minus)

    for _ in range(num_walks):
        position = np.array(start_pos, dtype=int)
        
        for step_count in range(max_steps):
            # Check if the walk hit the set A
            if tuple(position) in A:
                hits += 1
                break
                
            # Check if the walk escaped. np.linalg.norm is slow, so we check squared distance.
            if np.sum(position**2) > escape_radius**2:
                break

            # Take a random step
            delta = steps[np.random.randint(0, 2 * d)]
            position += delta
        # If the loop finishes without hitting or escaping, it's counted as an escape
        
    return hits / num_walks

def main():
    """
    Main function to run the simulation and demonstrate the property.
    """
    d = 3
    # A finite set like the origin is a classic example of a transient set in d>=3.
    A = {(0, 0, 0)}
    
    # Simulation parameters
    num_walks_per_point = 10000
    max_steps_per_walk = 5000
    # Escape radius should be significantly larger than the starting distances
    escape_radius = 50 

    print("The question asks if a set A can be transient if it's hit with probability 1 from infinitely many starting points.")
    print("The answer is No. This script demonstrates the reason: for any transient set, the hitting probability must approach 0 as the start point 'x' moves to infinity.")
    print("An infinite set of points where the hitting probability is 1 would contradict this.")
    print("\nThe following simulation illustrates this principle for a known transient set, A = {(0,0,0)} in 3D.")
    print("We estimate the hitting probability P_x(tau_A < infinity) from starting points 'x' at increasing distances from A.")
    print("-" * 55)
    print(f"Simulation Parameters:")
    print(f"  Dimension (d): {d}")
    print(f"  Transient Set (A): {A}")
    print(f"  Number of walks per point: {num_walks_per_point}")
    print(f"  Max steps per walk: {max_steps_per_walk}")
    print(f"  Escape radius: {escape_radius}")
    print("-" * 55)
    print("Distance | Start Point x     | Estimated P_x(hit A)")
    print("-" * 55)

    # Test for starting points at various distances from the origin
    for distance in range(1, 16):
        start_pos = tuple([distance] + [0] * (d - 1))
        
        prob = simulate_hitting_probability(
            d=d,
            A=A,
            start_pos=start_pos,
            num_walks=num_walks_per_point,
            max_steps=max_steps_per_walk,
            escape_radius=escape_radius
        )
            
        print(f"{distance:<8} | {str(start_pos):<17} | {prob:.4f}")
        
    print("-" * 55)
    print("\nAs the distance of the starting point 'x' increases, the estimated hitting probability")
    print("clearly decreases, approaching 0. This supports the theoretical argument.")
    print("Therefore, a set A with the property P_x(tau_A < infinity)=1 for infinitely many x cannot be transient.")

if __name__ == "__main__":
    main()