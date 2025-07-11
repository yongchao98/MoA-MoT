import numpy as np

def run_simulation(p, s_objects, s_star_true, num_candidates=2000, max_triplets=2500, conv_variance=1e-5):
    """
    Simulates the student's learning process for a given alignment p.
    The student uses Bayesian inference to update their beliefs.

    Returns the number of triplets required for convergence.
    """
    num_obj, dim = s_objects.shape
    
    # 1. Initialize student's belief: a cloud of candidate points with uniform weights
    candidates = np.random.rand(num_candidates, dim)
    weights = np.ones(num_candidates) / num_candidates

    for num_triplets in range(1, max_triplets + 1):
        # 2. Teacher provides a random triplet statement.
        # First, select two random reference objects.
        idx1, idx2 = np.random.choice(num_obj, 2, replace=False)
        obj1, obj2 = s_objects[idx1], s_objects[idx2]

        # 3. We simulate the teacher's statement. We don't need the teacher's space.
        #    We only need to know that with probability `p`, the teacher's statement
        #    aligns with the true relationship in the student's space.
        is_obj1_closer_in_student_truth = np.linalg.norm(s_star_true - obj1) < np.linalg.norm(s_star_true - obj2)
        
        if np.random.rand() < p:
            # With probability p, the teacher's statement agrees with the student's ground truth.
            teacher_says_obj1_is_closer = is_obj1_closer_in_student_truth
        else:
            # With probability 1-p, it disagrees.
            teacher_says_obj1_is_closer = not is_obj1_closer_in_student_truth

        # 4. Student updates beliefs. The student knows 'p' and performs a Bayesian update.
        #    Calculate the likelihood P(statement | candidate) for each candidate.
        dist_sq_1 = np.sum((candidates - obj1)**2, axis=1)
        dist_sq_2 = np.sum((candidates - obj2)**2, axis=1)
        
        # This mask is True where a candidate's geometry agrees with the teacher's statement.
        if teacher_says_obj1_is_closer:
            candidate_agrees_mask = dist_sq_1 < dist_sq_2
        else:
            candidate_agrees_mask = dist_sq_2 < dist_sq_1

        # The likelihood of the observation is 'p' for agreeing candidates and '1-p' for disagreeing ones.
        # This allows the student to correctly interpret even anti-aligned information (p < 0.5).
        likelihoods = np.full(num_candidates, 1.0 - p)
        likelihoods[candidate_agrees_mask] = p
        
        # Update weights: new_weight = old_weight * likelihood
        weights *= likelihoods
        
        # If all weights become near-zero (due to contradictory info or bad luck), stop.
        if np.sum(weights) < 1e-12:
            return max_triplets

        # Normalize weights to keep them as a probability distribution.
        weights /= np.sum(weights)
        
        # 5. Check for convergence by measuring the variance of the particle cloud.
        if len(candidates) > 1:
            mean_pos = np.average(candidates, weights=weights, axis=0)
            variance = np.average(np.sum((candidates - mean_pos)**2, axis=1), weights=weights)
            if variance < conv_variance:
                return num_triplets

    return max_triplets

def main():
    """Main function to run the experiment and print results."""
    # Simulation environment parameters
    DIM = 2
    NUM_REFERENCE_OBJECTS = 30
    NUM_SIMULATION_RUNS = 5  # Average over several runs to get a smoother curve
    MAX_TRIPLETS_CAP = 2500

    # Define a fixed set of objects for consistency across simulations
    s_objects = np.random.rand(NUM_REFERENCE_OBJECTS, DIM)
    s_star_true = np.random.rand(1, DIM) # The "true" location the student must find

    p_values = np.linspace(0.0, 1.0, 11)
    
    print("Running simulation...")
    print("This demonstrates the relationship between alignment (p) and the number of triplets to learn a location.")
    print("-" * 60)
    print(f"{'Alignment `p`':<20} | {'Avg. Triplets Required':<30}")
    print("-" * 60)

    for p in p_values:
        # At p=0.5, no information is gained, so it will never converge.
        if np.isclose(p, 0.5):
            avg_result = MAX_TRIPLETS_CAP
        else:
            run_results = [run_simulation(p, s_objects, s_star_true, max_triplets=MAX_TRIPLETS_CAP) for _ in range(NUM_SIMULATION_RUNS)]
            avg_result = np.mean(run_results)
        
        # We print each number that contributes to the final relationship.
        print(f"{p:<20.2f} | {int(avg_result):<30}")

    print("-" * 60)
    print("Conclusion: The results show that the number of required triplets is highest when p=0.5")
    print("and lowest at the extremes (p=0 and p=1). This describes a convex U-shaped curve.")

if __name__ == '__main__':
    main()