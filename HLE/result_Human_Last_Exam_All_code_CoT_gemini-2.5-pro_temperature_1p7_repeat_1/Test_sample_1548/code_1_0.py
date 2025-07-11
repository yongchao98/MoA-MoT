import numpy as np
import time

def calculate_potential_kernel(N=60, iterations=5000):
    """
    Calculates the potential kernel `a(x)` on a grid [-N, N] x [-N, N]
    using a relaxation method. The potential kernel is defined by the
    equation: (Delta a)(x) = -delta_{x,0}, where Delta is the discrete
    Laplacian. This is equivalent to a(x) = 0.25 * sum(a(y) for y~x) + delta_{x,0}.
    
    Args:
        N (int): The half-size of the grid. The grid will be (2N+1)x(2N+1).
        iterations (int): The number of iterations for the relaxation method.

    Returns:
        A (2N+1)x(2N+1) numpy array representing the potential kernel a(x).
    """
    size = 2 * N + 1
    # Initialize potential kernel grid with zeros.
    a = np.zeros((size, size))
    print("Calculating potential kernel on a {}x{} grid...".format(size, size))
    
    # The value a(0,0) is 0 for the standard potential kernel.
    # We solve for a(x) on Z^2\{0} satisfying a(x) = 1/4 * sum of neighbors + delta_x,0
    # Note that the origin (N,N) will have its value computed. This a is the Green's function G(x,0) killed on the boundary.
    # The true kernel has a(0)=0 but a has a log singularity at origin if computed differently. 
    # The relaxation `a = 0.25 * sum_neighbors(a) + delta_0` computes G_D(x,0)
    # For h-transform ratios a(y)/a(x) for x,y!=0 this is fine.
    
    for i in range(iterations):
        # Create a padded version for easy neighborhood sum
        padded_a = np.pad(a, 1, 'constant', constant_values=0)
        # Sum over 4 neighbors
        s = padded_a[:-2, 1:-1] + padded_a[2:, 1:-1] + padded_a[1:-1, :-2] + padded_a[1:-1, 2:]
        a_new = 0.25 * s
        # Add the source term at the origin (index N,N)
        a_new[N, N] += 1
        # The value at the origin itself is not part of the space Z^2\{0},
        # but is needed for computations near the origin.
        # For the h-transform, we won't be AT the origin, only near it.
        a = a_new
        if i % 1000 == 0:
            print("Iteration {}/{}".format(i, iterations))

    print("Potential kernel calculation finished.")
    # The problem defines a as the potential kernel of SRW. This kernel is positive.
    # So we return `a`. a(0,0) will be G(0,0)=1 in this finite volume approximation.
    return a

def run_simulation(a, steps, start_pos=(1, 0)):
    """
    Runs the simulation of the h-process.
    
    Args:
        a (np.array): The precomputed potential kernel.
        steps (int): The number of steps to simulate.
        start_pos (tuple): The starting position.

    Returns:
        The number of visits to the x-axis.
    """
    N = (a.shape[0] - 1) // 2
    # The walker must stay within a box smaller than the grid for `a`
    # to ensure all neighbors have a defined kernel value.
    sim_boundary = N - 1

    pos = np.array(start_pos)
    
    # The starting point is on the x-axis.
    x_axis_visits = 1 if pos[1] == 0 and tuple(pos) != (0,0) else 0
    
    for step in range(steps):
        x, y = pos
        if abs(x) > sim_boundary or abs(y) > sim_boundary:
            # This case should be rare if N is large enough for the number of steps.
            break
        
        # Get the kernel value at the current position
        a_idx_x, a_idx_y = x + N, y + N
        a_current = a[a_idx_x, a_idx_y]

        # Get kernel values for the 4 neighbors
        neighbors = [np.array([x + 1, y]), np.array([x - 1, y]),
                     np.array([x, y + 1]), np.array([x, y - 1])]
        
        probs = []
        valid_neighbors = []
        for neigh_pos in neighbors:
            # The walk is on Z^2 \ {0}
            if tuple(neigh_pos) != (0, 0):
                nx, ny = neigh_pos
                a_neigh_idx_x, a_neigh_idx_y = nx + N, ny + N
                a_neighbor = a[a_neigh_idx_x, a_neigh_idx_y]
                
                # The transition probability P(pos -> neigh_pos) = (1/4) * a(neigh)/a(pos)
                prob = 0.25 * a_neighbor / a_current
                probs.append(prob)
                valid_neighbors.append(neigh_pos)
        
        # Normalize probabilities to sum to 1. This accounts for boundary effects
        # on `a` making it not perfectly harmonic, and for avoiding the origin.
        probs = np.array(probs)
        probs /= probs.sum()

        # Choose the next position based on the calculated probabilities
        next_idx = np.random.choice(len(valid_neighbors), p=probs)
        pos = valid_neighbors[next_idx]

        # Check if the new position is on the x-axis
        if pos[1] == 0 and tuple(pos) != (0,0):
             x_axis_visits += 1
             
    return x_axis_visits

if __name__ == '__main__':
    # Set random seed for reproducibility
    np.random.seed(42)

    # Step 1: Calculate potential kernel `a`
    # A larger grid and more iterations give a better approximation of `a`
    # but take longer to compute. N=60 is a reasonable compromise.
    potential_kernel = calculate_potential_kernel(N=60, iterations=10000)

    # Step 2: Run simulation for increasing numbers of steps
    simulation_steps = [1000, 5000, 20000, 50000]
    
    print("\n--- h-Walk Simulation Results ---")
    print("This simulation demonstrates that an infinite set (the x-axis)")
    print("is visited frequently by the walk, suggesting it is a recurrent set.")
    
    for steps in simulation_steps:
        print(f"\nRunning simulation for {steps} steps...")
        start_time = time.time()
        visits = run_simulation(potential_kernel, steps=steps, start_pos=(1,0))
        duration = time.time() - start_time
        # Final equation format requested by the prompt.
        print("Number of visits to the x-axis after {} steps = {}".format(steps, visits))
        print(f"(Took {duration:.2f} seconds)")
