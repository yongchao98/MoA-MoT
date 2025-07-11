import numpy as np

def simulate_learning(p, confidence_threshold=10, max_triplets=50000):
    """
    Simulates how many triplets are needed for a student to make a confident decision
    based on a teacher's statements, given an alignment probability 'p'.

    Args:
        p (float): The probability of the teacher and student agreeing on a triplet.
        confidence_threshold (int): The vote difference needed for the student to be confident.
        max_triplets (int): A cutoff to prevent infinite loops when p is near 0.5.

    Returns:
        int: The number of triplets required to reach the confidence threshold.
    """
    # At p=0.5, there is no information, so it would take infinite triplets.
    # We return the maximum number to represent this for the simulation.
    if p == 0.5:
        return max_triplets

    # We can assume the ground truth is that o* is closer to o_j than o_k.
    # The teacher's statement will agree with this ground truth with probability p.
    prob_teacher_agrees_with_truth = p
    
    votes_for_truth = 0
    votes_against_truth = 0
    num_triplets = 0

    # The student collects triplets until one option has a clear lead.
    while abs(votes_for_truth - votes_against_truth) < confidence_threshold:
        if num_triplets >= max_triplets:
            # Reached max, this happens when p is very close to 0.5
            return max_triplets
            
        num_triplets += 1
        
        # Simulate a single triplet statement from the teacher
        if np.random.rand() < prob_teacher_agrees_with_truth:
            # Teacher's statement matches the truth
            votes_for_truth += 1
        else:
            # Teacher's statement contradicts the truth
            votes_against_truth += 1
            
    return num_triplets

def run_simulation():
    """
    Runs the simulation for a range of p values and prints the results.
    """
    # Define the alignment probabilities to test
    p_values = np.linspace(0, 1, 21)
    
    # Run multiple trials per p-value for a smoother average
    num_trials = 200
    
    # --- Output Results ---
    print("This script simulates the relationship between representational alignment (p)")
    print("and the number of triplets a teacher must send to a student.")
    print("\nA lower 'Avg. Triplets Needed' means learning is more efficient.\n")
    print("-" * 55)
    print(f"{'Alignment `p`':<20} | {'Avg. Triplets Needed':<25}")
    print("-" * 55)

    for p in p_values:
        # For p=0.5, the result is theoretically infinite.
        if p == 0.5:
            avg_triplets_str = "Infinity"
        else:
            total_triplets = 0
            for _ in range(num_trials):
                total_triplets += simulate_learning(p)
            avg_triplets = total_triplets / num_trials
            avg_triplets_str = f"{avg_triplets:.1f}"
            
        # Each line of this output shows a number in the relationship
        print(f"{p:<20.2f} | {avg_triplets_str:<25}")

    print("-" * 55)
    print("\nConclusion:")
    print("The number of triplets is minimal at the extremes (p=0 and p=1) and maximal")
    print("in the middle (p=0.5). This demonstrates a convex U-shaped relationship.")

if __name__ == "__main__":
    run_simulation()