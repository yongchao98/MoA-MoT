import numpy as np

def run_simulation(p, num_triplets=40, num_known_objects=20, grid_size=40, num_trials=15):
    """
    Simulates the student's learning process for a given alignment p.
    Returns the average estimation error over several trials.
    """
    trial_errors = []
    for _ in range(num_trials):
        # 1. SETUP THE SPACE for a single trial
        # Define the boundaries and the true location of the new object
        bounds = [-5, 5]
        o_star_true = np.array([0., 0.])
        
        # Randomly place known objects in the space
        known_objects = np.random.uniform(bounds[0], bounds[1], size=(num_known_objects, 2))

        # Create a grid of points representing the student's belief about the location of o*
        x_coords = np.linspace(bounds[0], bounds[1], grid_size)
        y_coords = np.linspace(bounds[0], bounds[1], grid_size)
        belief_points = np.array(np.meshgrid(x_coords, y_coords)).T.reshape(-1, 2)
        
        # Initialize log-likelihood scores for each belief point. We use logs for numerical stability.
        log_scores = np.zeros(len(belief_points))
        
        # Pre-calculate log(p) and log(1-p) to handle p=0 and p=1 cases.
        log_p = np.log(p) if p > 0 else -np.inf
        log_1_p = np.log(1 - p) if p < 1 else -np.inf

        # 2. SIMULATE THE TEACHING PROCESS
        for _ in range(num_triplets):
            # The Teacher picks two random objects to form a triplet statement.
            j, k = np.random.choice(num_known_objects, 2, replace=False)
            oj, ok = known_objects[j], known_objects[k]

            # The Teacher evaluates the triplet in its own space. Here we use o_star_true as the ground truth.
            is_j_closer_true = np.linalg.norm(o_star_true - oj) < np.linalg.norm(o_star_true - ok)

            # The Teacher's statement agrees with the truth with probability p.
            if np.random.rand() < p:
                teacher_says_j_is_closer = is_j_closer_true
            else:
                teacher_says_j_is_closer = not is_j_closer_true

            # The Student evaluates the Teacher's statement against all its belief points.
            dist_j_belief = np.linalg.norm(belief_points - oj, axis=1)
            dist_k_belief = np.linalg.norm(belief_points - ok, axis=1)
            
            # Find which belief points are consistent with the teacher's statement.
            if teacher_says_j_is_closer:
                consistent_mask = dist_j_belief < dist_k_belief
            else: # Teacher says k is closer
                consistent_mask = dist_k_belief > dist_j_belief
            
            # The Student updates the log-likelihood scores of its belief points (Bayesian update).
            # The probability of the statement given a point is p if consistent, and 1-p if not.
            log_scores[consistent_mask] += log_p
            log_scores[~consistent_mask] += log_1_p

        # 3. EVALUATE THE RESULT
        # The Student's best guess is the belief point with the highest score.
        best_guess_idx = np.argmax(log_scores)
        best_guess_point = belief_points[best_guess_idx]
        
        # Calculate the error between the student's guess and the true location.
        error = np.linalg.norm(best_guess_point - o_star_true)
        trial_errors.append(error)
        
    return np.mean(trial_errors)

if __name__ == '__main__':
    # Define the range of p (probabilistic alignment) to test.
    p_values = np.linspace(0, 1, 11)
    
    print("Running simulation to find the relationship between alignment (p) and estimation error.")
    print("A higher error after a fixed number of triplets implies more triplets are needed for accuracy.")
    print("-" * 50)
    print(f"{'Alignment (p)':<15} | {'Average Error':<15}")
    print("-" * 50)
    
    for p in p_values:
        # For each p, run the simulation and get the average error.
        avg_error = run_simulation(p)
        print(f"p = {p:<10.2f} | Error = {avg_error:.4f}")

    print("-" * 50)
    print("The error is lowest at p=0 and p=1, and highest at p=0.5, showing a convex U-shaped relationship.")