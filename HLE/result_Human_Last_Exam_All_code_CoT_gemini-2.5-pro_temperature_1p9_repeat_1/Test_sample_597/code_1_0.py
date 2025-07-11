import math
import random

def pcp_oracle(formula_type: str, proof: list[int]) -> float:
    """
    This is a hypothetical oracle representing the rejection probability of the Red/Blue PCP.
    Its behavior is what the problem statement allows us to assume. For this simulation, we
    need a placeholder for the "correct proof set". This function is not practically
    implementable as it would require solving the NP-complete problem itself.
    """
    proof_len = len(proof)
    if formula_type == "SATISFIABLE":
        # Let's imagine the single correct proof for a satisfiable formula is an alternating
        # sequence of 1s and 0s. A real proof structure would be far more complex.
        correct_proof = [(i + 1) % 2 for i in range(proof_len)]
        
        # Calculate relative Hamming distance from the correct proof.
        dist = sum(1 for i in range(proof_len) if proof[i] != correct_proof[i])
        relative_dist = dist / proof_len
        
        # The Red/Blue property states P_rej = Theta(delta). We'll use P_rej = 0.4 * delta.
        return 0.4 * relative_dist
    else:  # UNSATISFIABLE formula
        # The set of correct proofs is empty, so delta is 1 for any proof.
        # The Red property implies rejection probability must be a constant > 0.
        return 0.4

def estimate_rejection_prob(formula_type: str, proof: list[int]) -> float:
    """
    Estimates the rejection probability by querying the hypothetical oracle. In a real
    scenario, this would involve sampling the verifier a polynomial number of times.
    For our demonstration, we just call the idealized oracle.
    """
    return pcp_oracle(formula_type, proof)

def solve_np_with_red_blue_pcp(formula_type: str, proof_len: int):
    """
    Attempts to solve an NP-complete problem in randomized polynomial time by performing
    a local search, a method enabled by the existence of a Red/Blue PCP.
    """
    N = proof_len
    
    # 1. Start with an arbitrary proof string.
    current_proof = [0] * N
    current_rejection_prob = estimate_rejection_prob(formula_type, current_proof)

    print("Starting local search algorithm...")
    print(f"Initial estimated rejection prob: {current_rejection_prob:.6f}")

    # The loop will execute at most N times since the distance decreases with each pass.
    for pass_num in range(N):
        improved_in_pass = False
        best_neighbor = list(current_proof)
        best_neighbor_prob = current_rejection_prob
        
        # 2. Iterate through all N possible single-bit flips (the "neighborhood").
        for i in range(N):
            neighbor_proof = list(current_proof)
            neighbor_proof[i] = 1 - neighbor_proof[i]
            
            # 3. Estimate its rejection probability to see if it's an improvement.
            neighbor_prob = estimate_rejection_prob(formula_type, neighbor_proof)

            if neighbor_prob < best_neighbor_prob:
                best_neighbor = neighbor_proof
                best_neighbor_prob = neighbor_prob
                improved_in_pass = True

        # 4. If we found a better proof, move to it.
        if improved_in_pass:
            current_proof = best_neighbor
            current_rejection_prob = best_neighbor_prob
            print(f"Pass {pass_num + 1}: Found improvement. New estimated rejection prob: {current_rejection_prob:.6f}")
        else:
            # Otherwise, we have reached a minimum.
            print("No further improvement found. Local search terminated.")
            break
            
    # 5. Conclude based on the final rejection probability.
    print(f"\nFinal estimated rejection probability: {current_rejection_prob:.6f}")
    
    # We use a small threshold to check if the probability is effectively zero.
    if current_rejection_prob < 1.0 / (2 * N):
        print("Conclusion: Found a proof with (near) zero rejection probability. The formula is SATISFIABLE.")
    else:
        print("Conclusion: Final proof still has high rejection probability. The formula is UNSATISFIABLE.")


# --- Main execution showcasing the algorithm ---
PROOF_LENGTH = 64

print("--- Test Case 1: SATISFIABLE Formula ---")
solve_np_with_red_blue_pcp("SATISFIABLE", PROOF_LENGTH)

print("\n" + "="*50 + "\n")

print("--- Test Case 2: UNSATISFIABLE Formula ---")
solve_np_with_red_blue_pcp("UNSATISFIABLE", PROOF_LENGTH)