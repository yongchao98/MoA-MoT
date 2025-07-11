import math

def demonstrate_strategy_failure(M, prob_threshold=0.9):
    """
    Demonstrates that the proposed strategy fails to guarantee the desired success probability.

    Args:
        M (int): The number of boxes Alice randomizes her choice from.
        prob_threshold (float): The desired minimum probability of success.
    """
    print(f"Analyzing the strategy for M = {M} and a success probability threshold of {prob_threshold}.")

    # According to the strategy, the win probability P(Win) for a sequence 's' is:
    # P(Win|s) = 1 - |D_s intersect {1..M}| / M
    # For a guaranteed success rate of >= prob_threshold, we need:
    # |D_s intersect {1..M}| <= (1 - prob_threshold) * M for ALL sequences 's'.

    # Let's calculate the maximum allowed number of differing elements in the first M boxes.
    max_diffs_allowed = math.floor((1 - prob_threshold) * M)
    print(f"To guarantee success, the number of differences in the first {M} boxes must be <= {max_diffs_allowed}.")

    # An adversary can construct a sequence 's' that differs from its representative 'r'
    # on more than max_diffs_allowed positions.
    # Let's choose the number of differences to be one more than the allowed limit.
    adversarial_diffs = max_diffs_allowed + 1
    print(f"An adversary can construct a sequence that differs from its representative on the first {adversarial_diffs} positions.")

    # For such a sequence, the number of "losing" choices for Alice's unopened box K is adversarial_diffs.
    # The total number of choices is M.
    
    # Calculate the win probability for this adversarial sequence.
    win_prob = 1.0 - adversarial_diffs / M
    
    print("\nFor this adversarial sequence, the calculation is:")
    print(f"  Number of 'losing' choices for K = {adversarial_diffs}")
    print(f"  Total number of choices for K = {M}")
    print(f"  P(Win) = 1 - (Number of losing choices / Total choices)")
    print(f"  P(Win) = 1 - ({adversarial_diffs} / {M}) = {win_prob:.2f}")

    if win_prob < prob_threshold:
        print(f"\nThe calculated probability {win_prob:.2f} is less than the required threshold {prob_threshold}.")
        print("This demonstrates that the strategy is not guaranteed to succeed.")
    else:
        # This case should not be reached with this logic.
        print("\nThe strategy seems to work for this M, which contradicts the theoretical proof.")

if __name__ == '__main__':
    # Let's test with M=10 as in the problem description (9/10 probability)
    demonstrate_strategy_failure(M=10, prob_threshold=0.9)
    print("\n" + "="*50 + "\n")
    # Let's test with a larger M to show it doesn't help
    demonstrate_strategy_failure(M=100, prob_threshold=0.9)
