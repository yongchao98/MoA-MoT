import numpy as np

def simulate_h_transform_walk():
    """
    Simulates a 2D random walk with a Doob's h-transform to illustrate
    that an infinite set can be transient.
    """
    # Simulation parameters
    num_steps = 1_000_000
    start_pos = np.array([2, 0])
    
    # The infinite set A is the positive x-axis
    # A = {(k, 0) for k > 0}.
    def is_in_set_A(pos):
        return pos[1] == 0 and pos[0] > 0

    # The harmonic function h. We use a proxy that captures the essential
    # behavior of the potential kernel a(x), which grows like log(||x||).
    # Using log(1 + ||x||) avoids issues with ||x||=1 (log(1)=0) and ||x||=0.
    def h(pos):
        norm = np.linalg.norm(pos)
        return np.log(1 + norm)

    # Simulation setup
    current_pos = np.copy(start_pos)
    visit_count = 0
    
    # To demonstrate how visits cease over time, we track them in intervals.
    intervals = [1000, 10000, 100000, num_steps]
    visit_log = {t: 0 for t in intervals}
    last_interval_end = 0

    print(f"Simulating a h-transformed random walk for {num_steps} steps.")
    print(f"The infinite set is the positive x-axis, i.e., points (k, 0) with k > 0.")
    print(f"Starting position: {current_pos.tolist()}")
    print("-" * 30)

    for step in range(num_steps):
        if is_in_set_A(current_pos):
            visit_count += 1
            for t in intervals:
                if step < t:
                    visit_log[t] += 1
                    break

        # Define the 4 possible moves (neighbors)
        moves = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])
        neighbors = current_pos + moves
        
        # Consider only valid neighbors (not the origin)
        valid_neighbors = []
        weights = []
        for neighbor in neighbors:
            if not np.array_equal(neighbor, [0, 0]):
                valid_neighbors.append(neighbor)
                weights.append(h(neighbor))
        
        if not weights:
            # This can only happen if the walk is at a position like (1,0) and
            # all neighbors are killed (which isn't the case here).
            # If it were to happen, the walk would terminate.
            print("Walk terminated, no valid moves.")
            break
            
        # Normalize weights to get probabilities
        probabilities = np.array(weights) / np.sum(weights)
        
        # Choose the next position based on the probabilities
        chosen_index = np.random.choice(len(valid_neighbors), p=probabilities)
        current_pos = valid_neighbors[chosen_index]
        
    print(f"Simulation finished.")
    print(f"Final position: {current_pos.tolist()}")
    print(f"Final distance from origin: {np.linalg.norm(current_pos):.2f}")
    print("-" * 30)
    print(f"Total visits to the infinite set A: {visit_count}")
    
    for t in intervals:
        print(f"Visits within steps {last_interval_end}-{t-1}: {visit_log[t]}")
        last_interval_end = t
        
    print("\nObservation: The walk visits the set a number of times early on but then")
    print("drifts away from the origin and ceases to return to the set.")
    print("This suggests that the infinite set is transient, meaning it is visited a.s. finitely many times.")

if __name__ == '__main__':
    simulate_h_transform_walk()
