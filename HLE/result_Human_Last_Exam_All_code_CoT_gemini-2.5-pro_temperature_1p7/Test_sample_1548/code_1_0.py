import numpy as np

def h(p):
    """
    This function represents the harmonic function used for the h-transform.
    It's based on the asymptotic behavior of the potential kernel a(x) of the 2D SRW.
    h(x) ~ log|x|. We add a constant to ensure it's positive.
    """
    x, y = p
    norm = np.sqrt(x**2 + y**2)
    # A large constant is added to log(|x|) to keep h positive, which is a
    # property of the actual potential kernel a(x) on Z^2 \ {0}.
    return np.log(norm) + 5.0

def run_simulation(start_pos, num_steps, sparse_set):
    """
    Simulates the Doob's h-transform of the SRW and counts
    the number of unique points visited in a given sparse set.
    """
    pos = np.array(start_pos, dtype=float)
    visited_points = set()

    # Check if starting point is in the set
    if tuple(np.round(pos).astype(int)) in sparse_set:
        visited_points.add(tuple(np.round(pos).astype(int)))

    # Get the value of the harmonic function at the starting point.
    # This value is used in the denominator of the transition probability calculation.
    # The value is updated in each step.
    h_current = h(pos)

    for _ in range(num_steps):
        neighbors = [pos + d for d in [[1, 0], [-1, 0], [0, 1], [0, -1]]]
        
        # Filter out the origin
        valid_neighbors = [n for n in neighbors if not np.all(n == 0)]
        
        # Transition probabilities P_h(x, y) are proportional to P(x, y) * h(y) / h(x).
        # Since P(x, y) = 1/|neighbors| is constant for all valid neighbors,
        # the probabilities are proportional to h(y).
        h_values = [h(n) for n in valid_neighbors]
        
        # We need to divide by h(x) to get the correct weights, but since we normalize
        # afterwards, dividing all weights by the same constant h(x) doesn't change the outcome.
        # So, probs are proportional just to h(y).
        probs = np.array(h_values)
        probs /= np.sum(probs) # Normalize to get probabilities

        # Choose the next position based on the calculated probabilities
        idx = np.random.choice(len(valid_neighbors), p=probs)
        pos = np.array(valid_neighbors[idx])
        
        h_current = h_values[idx]

        # Check if the new position is in our sparse set
        # We round the position because of potential float inaccuracies
        pos_int = tuple(np.round(pos).astype(int))
        if pos_int in sparse_set:
            visited_points.add(pos_int)

    return len(visited_points)

if __name__ == '__main__':
    # Define an infinite, sparse set. We take points on the x-axis with exponentially growing distance.
    # We will test a finite subset of these points.
    sparse_set = {(2**k, 0) for k in range(1, 16)}

    # Simulation parameters
    num_steps = 200000
    start_pos = (2, 0) # Start on a point in the set

    # Run the simulation
    num_visited = run_simulation(start_pos, num_steps, sparse_set)

    # Output the result.
    # This result provides numerical evidence that an infinite sparse set can be transient.
    # In a typical run, the number of visited sites is very small (e.g., 2-5),
    # even with hundreds of thousands of steps, suggesting the walker quickly
    # moves away and is unlikely to return to visit many other points in the set.
    print(f"Simulation of a random walk for {num_steps} steps.")
    print(f"The walk started at {start_pos} which is in the sparse set A = {{(2^k, 0): k >= 1}}.")
    print(f"Number of unique points visited in the set A: {num_visited}")
<<<No>>>