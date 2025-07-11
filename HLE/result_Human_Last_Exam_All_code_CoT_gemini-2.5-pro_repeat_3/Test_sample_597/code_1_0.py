import random

def illustrate_blue_pcp_property():
    """
    This function simulates a toy PCP verifier to illustrate the Blue PCP property.
    The Blue property states that the rejection probability P_rej is at most
    proportional to the proof's distance from correctness (delta), i.e., P_rej <= C * delta.
    For a verifier with query complexity q, C is typically q.
    """
    # --- Parameters ---
    # The length of the proof string.
    proof_length = 20000
    # The verifier's query complexity (a small constant).
    q = 5
    # The relative Hamming distance of the provided proof from a correct one.
    delta = 0.1

    print("This script illustrates the 'Blue PCP' property: P_rej <= C * delta")
    print("----------------------------------------------------------------\n")
    print(f"Setting up a toy scenario:")
    print(f"Proof length N = {proof_length}")
    print(f"Verifier query complexity q = {q}. We will use C = q.")

    # --- Proof Setup ---
    # A "correct" proof (e.g., all zeros).
    pi_correct = [0] * proof_length
    # The number of errors in the provided proof.
    num_errors = int(delta * proof_length)
    print(f"Proof's relative distance from correct, delta = {delta} ({num_errors} errors)")

    # Create the provided proof 'pi' by introducing 'num_errors' into the correct proof.
    pi = list(pi_correct)
    error_indices = random.sample(range(proof_length), num_errors)
    for i in error_indices:
        pi[i] = 1 # Introduce an error.

    # --- Verifier Simulation ---
    num_simulations = 50000
    rejection_count = 0
    print(f"\nSimulating the verifier {num_simulations} times...")

    for _ in range(num_simulations):
        # The verifier randomly selects q positions to query.
        query_indices = random.sample(range(proof_length), q)

        # Our toy verifier rejects if it finds any error (a '1') in the queried positions.
        found_error = False
        for i in query_indices:
            if pi[i] == 1:
                found_error = True
                break
        if found_error:
            rejection_count += 1

    # --- Results ---
    # The observed rejection probability.
    observed_p_rej = rejection_count / num_simulations
    # The theoretical upper bound from the Blue property.
    theoretical_bound = q * delta

    print("\n--- Results ---")
    print(f"Observed rejection probability (P_rej): {observed_p_rej:.4f}")
    
    # The final equation is the check of the Blue property.
    print("\nFinal Equation Check (P_rej <= q * delta):")
    # We print each number in the equation.
    print(f"{observed_p_rej:.4f} <= {q} * {delta} = {theoretical_bound:.4f}")
    is_satisfied = observed_p_rej <= theoretical_bound
    print(f"\nIs the Blue property satisfied in this simulation? {'Yes.' if is_satisfied else 'No.'}")

if __name__ == "__main__":
    illustrate_blue_pcp_property()