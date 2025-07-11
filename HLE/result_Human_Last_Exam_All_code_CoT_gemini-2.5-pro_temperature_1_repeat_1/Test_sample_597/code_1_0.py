import math
import random

# --- Helper functions ---
def relative_hamming_distance(str1, str2):
    """Calculates the relative Hamming distance between two strings."""
    if len(str1) != len(str2):
        raise ValueError("Strings must have the same length.")
    if len(str1) == 0:
        return 0.0
    distance = sum(c1 != c2 for c1, c2 in zip(str1, str2))
    return distance / len(str1)

def distance_to_set(pi, pi_set):
    """Calculates the relative Hamming distance of a string pi to a set of strings pi_set."""
    if not pi_set:  # If the set is empty (NO instance)
        return 1.0
    return min(relative_hamming_distance(pi, correct_pi) for correct_pi in pi_set)

# --- Simulation of the Hypothetical PCP ---

class HypotheticalRedBluePCPVerifier:
    """
    A simulated verifier for a hypothetical PCP that is both Red and Blue.
    Its rejection probability is directly proportional to the distance from a correct proof.
    """
    def __init__(self, x, np_problem_checker):
        self.x = x
        self.pi_set = np_problem_checker(x)
        # The constants for the Red/Blue property: c*delta <= P[rej] <= C*delta
        # For simplicity, we simulate the ideal case where P[rej] = k * delta.
        self.PROPORTIONALITY_CONSTANT = 0.5 # An arbitrary constant between 0 and 1.

    def get_rejection_prob(self, pi):
        """
        Calculates the rejection probability for a given proof 'pi'.
        This simulates the core property of the Red/Blue PCP.
        """
        delta = distance_to_set(pi, self.pi_set)
        # The rejection probability is Theta(delta). We simulate it as k * delta.
        # P[reject] = k * delta(pi, Pi(x))
        rejection_prob = self.PROPORTIONALITY_CONSTANT * delta
        return rejection_prob

    def query(self, pi):
        """Simulates a single run of the verifier."""
        prob = self.get_rejection_prob(pi)
        return random.random() < prob

# --- A Dummy NP Problem for Demonstration ---
# For our demo, we use a problem in P, but pretend it's NP-complete.
# Problem: Does a binary string have exactly one '1'?
# A "correct proof" for a YES instance 'x' is the string 'x' itself.
def single_one_checker(x):
    """
    Checks if the input string 'x' is a YES or NO instance.
    Returns the set of correct proofs Pi(x).
    """
    if x.count('1') == 1:
        return {x}  # YES instance, Pi(x) = {x}
    else:
        return set()  # NO instance, Pi(x) is empty

# --- The Polynomial-Time Solver for NP ---

def solve_np_in_p_with_hypothetical_pcp(x, problem_checker, proof_len):
    """
    This function demonstrates how to solve an NP problem in polynomial time
    if a Red/Blue PCP for it existed.
    """
    print(f"--- Attempting to solve for input x = '{x}' ---")
    verifier = HypotheticalRedBluePCPVerifier(x, problem_checker)

    # To estimate probability, we run the verifier many times.
    # The number of samples needs to be polynomial in the proof length.
    # A change of 1 bit changes distance by 1/proof_len. We need to detect this.
    # O((proof_len)^2) samples are sufficient.
    num_samples = proof_len * proof_len * 10

    def estimate_rejection_prob(pi):
        """Helper to estimate rejection probability by sampling."""
        rejections = sum(1 for _ in range(num_samples) if verifier.query(pi))
        return rejections / num_samples

    # Start with a random or zero-proof.
    current_pi_list = ['0'] * proof_len
    current_pi = "".join(current_pi_list)
    
    # Greedily "self-correct" the proof to find a valid one if it exists.
    # This works because the verifier's rejection probability is a reliable
    # proxy for the actual distance to a correct proof.
    print("Starting greedy self-correction of the proof...")
    # We iterate enough times to ensure convergence. proof_len passes are sufficient.
    for pass_num in range(proof_len):
        made_a_change = False
        # Get the baseline rejection probability for the current proof.
        base_prob = estimate_rejection_prob(current_pi)
        
        for i in range(proof_len):
            # Try flipping a bit
            original_char = current_pi_list[i]
            current_pi_list[i] = '1' if original_char == '0' else '0'
            flipped_pi = "".join(current_pi_list)

            # Check if flipping the bit improved the proof (reduced rejection prob)
            flipped_prob = estimate_rejection_prob(flipped_pi)

            if flipped_prob < base_prob:
                # The flip was an improvement, so we keep it.
                current_pi = flipped_pi
                base_prob = flipped_prob
                made_a_change = True
                print(f"  Pass {pass_num+1}, Bit {i}: Flipped bit. New proof: {current_pi}, Est. Rej. Prob: {base_prob:.4f}")
            else:
                # Revert the flip if it wasn't an improvement.
                current_pi_list[i] = original_char
        
        if not made_a_change:
            print(f"  No improvement found in pass {pass_num+1}. Converged.")
            break
            
    final_pi = current_pi
    print(f"\nFinal candidate proof after self-correction: '{final_pi}'")

    # Final check: A correct proof has a rejection probability of 0.
    # A proof for a NO instance has a rejection probability > 0.
    final_rejection_prob = verifier.get_rejection_prob(final_pi)
    
    print(f"Final exact rejection probability for candidate proof is P[reject] = {verifier.PROPORTIONALITY_CONSTANT} * delta")
    print(f"  delta = distance('{final_pi}', Pi('{x}')) = {distance_to_set(final_pi, verifier.pi_set)}")
    print(f"  P[reject] = {verifier.PROPORTIONALITY_CONSTANT} * {distance_to_set(final_pi, verifier.pi_set)} = {final_rejection_prob}")

    # The Decision
    if final_rejection_prob == 0:
        print("\nConclusion: The solver found a proof with zero rejection probability.")
        print(f"This implies '{x}' is a YES instance.")
        return "YES"
    else:
        print("\nConclusion: The solver could not find a proof with zero rejection probability.")
        print(f"This implies '{x}' is a NO instance.")
        return "NO"

if __name__ == '__main__':
    # --- Run the simulation ---
    # Case 1: A YES instance for our dummy NP problem
    yes_instance = "00100"
    solve_np_in_p_with_hypothetical_pcp(yes_instance, single_one_checker, len(yes_instance))

    print("\n" + "="*70 + "\n")

    # Case 2: A NO instance for our dummy NP problem
    no_instance = "01100"
    solve_np_in_p_with_hypothetical_pcp(no_instance, single_one_checker, len(no_instance))
    
    print("\n" + "="*70 + "\n")
    print("Final Analysis:")
    print("The simulation shows that if a PCP were both Red and Blue, its verifier could be used to create a greedy algorithm that finds a correct proof in polynomial time.")
    print("Such an algorithm would solve any problem in NP in polynomial time.")
    print("This would mean P = NP, which contradicts the initial assumption of the problem.")
    print("Therefore, it is not possible for NP to have a PCP that is both Red and Blue (assuming P != NP).")
    print("\n<<<No>>>")
