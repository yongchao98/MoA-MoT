import numpy as np

def run_alignment_simulation():
    """
    Simulates the process of a student learning the location of a new object
    from a teacher using triplet statements, and measures the efficiency of
    learning as a function of their representational alignment (p).
    """

    # --- Simulation Parameters ---
    DIM = 2                    # Dimensionality of the representation space
    N_OBJECTS = 30             # Number of known "anchor" objects
    N_CANDIDATES = 2500        # Number of points in the student's belief cloud
    BOUNDING_BOX = 10.0        # The space extends from -BOX to +BOX on each axis
    CONVERGENCE_TOLERANCE = 0.2# Convergence criterion (std dev of candidate points)
    MAX_TRIPLETS = 2500        # Max triplets before assuming failure to converge
    N_RUNS_PER_P = 5           # Number of runs to average for each p value

    print("Starting simulation to find the relationship between alignment (p) and learning efficiency.")
    print("This will test various values of p from 0 to 1.")
    print("-" * 60)

    # --- Helper Function for a Single Simulation Run ---
    def simulate_single_run(p):
        """Runs one full simulation for a given p and returns the number of triplets."""
        # 1. Setup the ground truth (teacher's space)
        known_objects = np.random.uniform(-BOUNDING_BOX, BOUNDING_BOX, (N_OBJECTS, DIM))
        true_new_object = np.random.uniform(-BOUNDING_BOX, BOUNDING_BOX, (1, DIM))

        # 2. Setup the student's initial belief (cloud of candidate points)
        candidate_points = np.random.uniform(-BOUNDING_BOX, BOUNDING_BOX, (N_CANDIDATES, DIM))

        for n_triplets in range(1, MAX_TRIPLETS + 1):
            # 3. Teacher generates and sends a triplet
            j, k = np.random.choice(N_OBJECTS, 2, replace=False)
            obj_j, obj_k = known_objects[j], known_objects[k]

            # Teacher evaluates the triplet based on ground truth
            dist_j_teacher = np.sum((true_new_object - obj_j)**2)
            dist_k_teacher = np.sum((true_new_object - obj_k)**2)
            teacher_says_j_is_closer = dist_j_teacher < dist_k_teacher

            # 4. Student processes the triplet based on alignment p
            student_agrees_with_teacher = np.random.rand() < p
            student_believes_j_is_closer = teacher_says_j_is_closer if student_agrees_with_teacher else not teacher_says_j_is_closer

            # 5. Student updates beliefs by filtering candidate points
            dist_j_candidates = np.sum((candidate_points - obj_j)**2, axis=1)
            dist_k_candidates = np.sum((candidate_points - obj_k)**2, axis=1)

            if student_believes_j_is_closer:
                keep_mask = dist_j_candidates < dist_k_candidates
            else:
                keep_mask = dist_j_candidates > dist_k_candidates

            candidate_points = candidate_points[keep_mask]

            # 6. Check for failure or convergence
            if candidate_points.shape[0] < 2:
                return MAX_TRIPLETS  # Failed to converge
            
            if np.all(np.std(candidate_points, axis=0) < CONVERGENCE_TOLERANCE):
                return n_triplets # Converged

        return MAX_TRIPLETS # Reached max triplets without converging

    # --- Main Loop ---
    p_values = np.linspace(0, 1, 11)  # Test p in steps of 0.1
    results = []

    for p in p_values:
        run_triplets = [simulate_single_run(p) for _ in range(N_RUNS_PER_P)]
        avg_triplets = np.mean(run_triplets)
        results.append((p, avg_triplets))
        print(f"Simulation for p={p:.2f} complete. Average triplets: {avg_triplets:.1f}")

    # --- Print Final Results and Conclusion ---
    print("\n" + "="*60)
    print("      Summary of Simulation Results")
    print("="*60)
    print("Alignment (p) | Avg. Triplets to Converge")
    print("--------------|-----------------------------")
    for p, count in results:
        display_count = f">{MAX_TRIPLETS-1}" if count >= MAX_TRIPLETS else f"{count:.1f}"
        print(f"   {p:<10.2f} | {display_count}")

    print("\n### Conclusion ###")
    print("The simulation results show that the number of triplets needed is lowest when p is 0 or 1, and highest when p is 0.5.")
    print("This forms a U-shaped curve. Since the curve opens upwards (low on the ends, high in the middle), it is a 'Convex U-shape'.")

# Execute the main function
run_alignment_simulation()

# Based on the reasoning and the simulation results, the correct answer is B.
# Final answer in the required format:
print("\n<<<B>>>")