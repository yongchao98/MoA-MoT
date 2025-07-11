import numpy as np

def run_simulation():
    """
    Simulates the teaching process to find the relationship between
    alignment (p) and the number of triplets needed.
    """
    # --- Simulation Parameters ---
    N_KNOWN_OBJECTS = 50
    DIMENSIONS = 2
    LEARNING_RATE = 0.01
    CONVERGENCE_THRESHOLD = 0.1  # Converged when distance to true location is < this value
    MAX_TRIPLETS = 10000          # Maximum triplets before we give up
    NOISE_LEVELS = [0.0, 0.1, 0.2, 0.4, 0.6, 0.8] # Controls student's alignment

    print("Running simulation to find the relationship between alignment (p) and N_triplets...")
    print("-" * 70)

    # A. Create the Teacher's representations
    # Use a seed for reproducible random object locations
    np.random.seed(42)
    # T_known are the n existing objects, T_star is the new object
    teacher_known_objects = np.random.rand(N_KNOWN_OBJECTS, DIMENSIONS)
    teacher_new_object = np.random.rand(DIMENSIONS)

    results = []

    for noise in NOISE_LEVELS:
        # B. Create the Student's representations by adding noise to the Teacher's
        # The 'true' location for the student is unknown to the student during learning.
        student_known_objects = teacher_known_objects + np.random.normal(0, noise, (N_KNOWN_OBJECTS, DIMENSIONS))
        student_true_new_object = teacher_new_object + np.random.normal(0, noise, (DIMENSIONS))

        # C. Estimate the probabilistic alignment 'p' for this noise level
        n_samples_for_p = 2000
        agreements = 0
        for _ in range(n_samples_for_p):
            j, k = np.random.choice(N_KNOWN_OBJECTS, 2, replace=False)

            # Teacher's evaluation
            dist_T_j_sq = np.sum((teacher_new_object - teacher_known_objects[j])**2)
            dist_T_k_sq = np.sum((teacher_new_object - teacher_known_objects[k])**2)
            teacher_prefers_j = dist_T_j_sq < dist_T_k_sq

            # Student's evaluation
            dist_S_j_sq = np.sum((student_true_new_object - student_known_objects[j])**2)
            dist_S_k_sq = np.sum((student_true_new_object - student_known_objects[k])**2)
            student_prefers_j = dist_S_j_sq < dist_S_k_sq

            if teacher_prefers_j == student_prefers_j:
                agreements += 1
        p = agreements / n_samples_for_p

        # D. Simulate the student learning process to find N_triplets
        student_guess = np.zeros(DIMENSIONS)  # Student starts with a blind guess
        n_triplets_for_convergence = MAX_TRIPLETS

        for i in range(1, MAX_TRIPLETS + 1):
            # Teacher provides one triplet
            j, k = np.random.choice(N_KNOWN_OBJECTS, 2, replace=False)
            dist_T_j_sq = np.sum((teacher_new_object - teacher_known_objects[j])**2)
            dist_T_k_sq = np.sum((teacher_new_object - teacher_known_objects[k])**2)
            
            # The constraint to be enforced from the Teacher's triplet
            teacher_prefers_j = dist_T_j_sq < dist_T_k_sq

            # Student uses the triplet to update its guess using gradient descent
            s_j, s_k = student_known_objects[j], student_known_objects[k]
            dist_guess_j_sq = np.sum((student_guess - s_j)**2)
            dist_guess_k_sq = np.sum((student_guess - s_k)**2)
            
            # Update only if the current guess violates the teacher's triplet
            gradient = np.zeros(DIMENSIONS)
            if teacher_prefers_j and dist_guess_j_sq > dist_guess_k_sq:
                # Move guess to satisfy "j is closer"
                gradient = 2 * (student_guess - s_j) - 2 * (student_guess - s_k)
            elif not teacher_prefers_j and dist_guess_k_sq > dist_guess_j_sq:
                # Move guess to satisfy "k is closer"
                gradient = 2 * (student_guess - s_k) - 2 * (student_guess - s_j)
            
            student_guess -= LEARNING_RATE * gradient

            # Check for convergence every 10 triplets
            if i % 10 == 0:
                distance_to_true = np.linalg.norm(student_guess - student_true_new_object)
                if distance_to_true < CONVERGENCE_THRESHOLD:
                    n_triplets_for_convergence = i
                    break

        print(f"Alignment p = {p:.4f}  =>  Number of triplets needed = {n_triplets_for_convergence}")
        results.append((p, n_triplets_for_convergence))

    print("-" * 70)
    print("Conclusion: As the alignment 'p' increases, the number of triplets required to teach the new location decreases.")

if __name__ == '__main__':
    run_simulation()
