import random

def get_representative_A(sequence):
    """
    For case (A), all sequences are in the same equivalence class.
    We can choose the all-zero sequence as the representative.
    This choice is fixed and known to the adversary.
    """
    return [0] * len(sequence)

def run_alice_strategy(sequence, K, get_representative_func):
    """
    Simulates Alice's strategy and calculates her success probability for a given sequence.

    Args:
        sequence: The list of numbers in the boxes.
        K: The number of boxes Alice leaves closed.
        get_representative_func: A function to get the representative for a sequence.

    Returns:
        The probability of success for Alice.
    """
    # Alice leaves the first K boxes closed and observes the rest.
    # The tail of the sequence determines the equivalence class and thus the representative.
    representative = get_representative_func(sequence)

    # Alice chooses one of the K closed boxes at random and guesses its value
    # to be the corresponding value from the representative sequence.
    # We calculate the number of boxes among the first K where her guess would be wrong.
    closed_boxes_indices = range(K)
    num_wrong_guesses = len({i for i in closed_boxes_indices if sequence[i] != representative[i]})

    # The probability of picking a "good" box is (K - num_wrong_guesses) / K.
    prob_success = (K - num_wrong_guesses) / K

    print(f"For the given sequence, Alice leaves the first {K} boxes closed.")
    print(f"Number of boxes in this set where her guess would be wrong: {num_wrong_guesses}")
    print(f"Her probability of success is ({K} - {num_wrong_guesses}) / {K} = {prob_success:.2f}")
    
    return prob_success

# We need a dummy representative function for Case (B) to simulate the non-constructive choice.
# The key property is that sequences with the same tail have the same representative.
REPRESENTATIVE_CACHE = {}
def get_representative_B(sequence):
    """
    A dummy representative function for Case (B).
    It ensures that sequences with the same tail have the same representative.
    """
    K = 100 # Alice's strategy parameter
    tail = tuple(sequence[K:])
    if tail not in REPRESENTATIVE_CACHE:
        # For simulation, we can just use the sequence itself as the first representative of its class.
        REPRESENTATIVE_CACHE[tail] = list(sequence)
    return REPRESENTATIVE_CACHE[tail]

def demonstrate_failure():
    """Demonstrates that the strategy fails for both cases."""
    K = 100
    required_prob = 0.9
    num_errors_to_force_failure = int(K * (1 - required_prob)) + 1
    
    print("--- Analysis of Alice's Strategy ---")
    print(f"Alice chooses to leave K={K} boxes closed.")
    print(f"To guarantee >= {required_prob*100}% success, her guess must be correct for at least {K - num_errors_to_force_failure +1} of the {K} boxes.")
    print(f"This means the number of 'error' boxes must be <= {num_errors_to_force_failure - 1}.")
    print(f"An adversary can win by forcing {num_errors_to_force_failure} or more errors.\n")

    # --- Case (A) ---
    print("--- Case (A): Sequences are eventually zero ---")
    adversarial_sequence_A = [1] * num_errors_to_force_failure + [0] * (K - num_errors_to_force_failure + 100)
    print(f"Adversary chooses a sequence with {num_errors_to_force_failure} non-zero values in the first {K} boxes.")
    prob_A = run_alice_strategy(adversarial_sequence_A, K, get_representative_A)
    if prob_A < required_prob:
        print("Conclusion: Alice's success probability is less than 9/10. The strategy fails for case (A).\n")
    
    # --- Case (B) ---
    print("--- Case (B): No assumptions on the sequence ---")
    base_sequence = [0] * (K + 100)
    s_star = get_representative_B(base_sequence) # Adversary doesn't know s_star, but knows it exists.
    
    adversarial_sequence_B = list(s_star) # Start with the representative
    # Modify it in num_errors_to_force_failure positions to create differences
    for i in range(num_errors_to_force_failure):
        adversarial_sequence_B[i] = s_star[i] + 1
    
    print(f"Adversary creates a sequence that differs from its representative at {num_errors_to_force_failure} positions in the first {K} boxes.")
    prob_B = run_alice_strategy(adversarial_sequence_B, K, get_representative_B)
    if prob_B < required_prob:
        print("Conclusion: Alice's success probability is less than 9/10. The strategy fails for case (B).")

demonstrate_failure()