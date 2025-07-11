import math

class HypotheticalRedBluePCP:
    """
    A mock class representing the hypothetical Red and Blue PCP verifier for 3-SAT.
    This class simulates the verifier's behavior to demonstrate the argument.
    It cannot be implemented in reality without solving P vs NP.
    """
    def __init__(self, formula, is_satisfiable):
        """
        Initializes the mock PCP system.
        - formula: A representation of the 3-SAT instance.
        - is_satisfiable: A flag to tell the mock verifier how to behave.
        """
        # For demonstration, the proof length is based on the number of variables.
        # In a real PCP, this would be a specific polynomial.
        num_vars = self._count_vars(formula)
        self.proof_length = num_vars * num_vars if num_vars > 0 else 4

        # If the formula is satisfiable, there is a non-empty set of correct proofs Π(x).
        # We simplify by assuming the all-ones proof is the *only* correct proof.
        # This is sufficient to demonstrate the core logic of the distance-based search.
        if is_satisfiable:
            # A dummy correct proof (all ones)
            self.correct_proofs = [[1] * self.proof_length]
        else:
            # If not satisfiable, the set of correct proofs is empty.
            self.correct_proofs = []

        # The constants for the Red and Blue property: c*δ <= P(reject) <= C*δ
        # We'll use the midpoint for our simulation.
        self.theta_constant = 0.5

    def _count_vars(self, formula):
        """Helper to get a measure of formula size."""
        if not formula: return 0
        all_vars = set()
        for clause in formula:
            for literal in clause:
                all_vars.add(abs(literal))
        return len(all_vars)

    def _relative_hamming_distance(self, proof):
        """Calculates δ(proof, Π(x))."""
        if not self.correct_proofs:
            # Per the problem, distance from the empty set is 1.
            return 1.0

        min_dist = float('inf')
        for correct_proof in self.correct_proofs:
            dist = sum(1 for i in range(self.proof_length) if proof[i] != correct_proof[i])
            min_dist = min(min_dist, dist)

        return min_dist / self.proof_length

    def get_rejection_probability(self, proof):
        """
        Calculates P(reject) = Θ(δ). Our mock uses a simple linear relation.
        """
        delta = self._relative_hamming_distance(proof)
        return self.theta_constant * delta

    def estimate_distance(self, proof):
        """
        Simulates the process of estimating the distance by sampling the verifier.
        A real algorithm would run the verifier many times and use the rejection
        rate to estimate P(reject), then divide by the Θ constant.
        Our simulation calculates it directly for simplicity.
        """
        p_reject = self.get_rejection_probability(proof)
        estimated_delta = p_reject / self.theta_constant
        return estimated_delta


def solve_sat_with_hypothetical_pcp(formula, is_satisfiable_for_mock):
    """
    Implements the polynomial-time 3-SAT solver using the hypothetical PCP.
    """
    print("--- Running Polynomial-Time SAT Solver ---")
    print(f"Formula: {formula}")

    pcp = HypotheticalRedBluePCP(formula, is_satisfiable_for_mock)
    m = pcp.proof_length
    
    # The theory guarantees convergence in at most `m` steps of improvement.
    max_iterations = m + 1 

    # Start with an arbitrary proof (e.g., all zeros)
    current_proof = [0] * m
    print(f"Proof length (m): {m}")
    print("Starting local search with an all-zeros proof.")

    for i in range(max_iterations):
        # Estimate the distance of the current proof
        current_dist = pcp.estimate_distance(current_proof)
        
        # This print statement shows the key number in our process: the distance.
        print(f"\nIteration {i+1}: Current estimated distance δ = {current_dist:.4f}")

        # If distance is effectively zero, we have found a correct proof.
        if current_dist < 1.0 / (2 * m):
            print("Distance is close to zero. Conclusion: Formula is SATISFIABLE.")
            return True

        # Perform greedy local search: find the single-bit flip that most reduces distance.
        best_flip_idx = -1
        best_flip_dist = current_dist

        for j in range(m):
            flipped_proof = list(current_proof)
            flipped_proof[j] = 1 - flipped_proof[j]
            
            dist = pcp.estimate_distance(flipped_proof)
            
            if dist < best_flip_dist:
                best_flip_dist = dist
                best_flip_idx = j

        # Update the proof if a better one was found
        if best_flip_idx != -1:
            print(f"Found an improvement. Flipping bit {best_flip_idx} to reduce distance.")
            current_proof[best_flip_idx] = 1 - current_proof[best_flip_idx]
        else:
            # Theory shows this only happens if δ > 0 for NO instances, or δ = 0 for YES instances.
            # Since we already checked for δ=0, being here implies a NO instance.
            print("No single-bit flip improves the distance. Conclusion: Formula is UNSATISFIABLE.")
            return False

    print("\nReached max iterations without finding a proof. Conclusion: Formula is UNSATISFIABLE.")
    return False

# --- Main execution to demonstrate the algorithm ---

# Case 1: A satisfiable formula
sat_formula = [[1, -2], [-1, 2]]
solve_sat_with_hypothetical_pcp(sat_formula, is_satisfiable_for_mock=True)

print("\n" + "="*50 + "\n")

# Case 2: An unsatisfiable formula
unsat_formula = [[1], [-1]]
solve_sat_with_hypothetical_pcp(unsat_formula, is_satisfiable_for_mock=False)
