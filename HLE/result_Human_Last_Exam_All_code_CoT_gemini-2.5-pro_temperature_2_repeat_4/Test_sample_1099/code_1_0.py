import math

def calculate_simulation_resources():
    """
    Calculates and prints the minimal resources needed to simulate the correlations
    of a singlet state with a Local Hidden Variable (LHV) model.
    """

    print("--- Calculating Minimal Resources for Simulating Singlet State Correlations ---")
    print("\nThis model calculates the resources needed to classically simulate the correlations")
    print("obtained from POVM measurements on a bipartite singlet quantum state.\n")

    # --- Resource 1: Classical Communication ---
    communication_cost = 1
    print("="*60)
    print("Resource 1: Classical Communication")
    print("="*60)
    print("To perfectly simulate the correlations for all projective measurements, the minimal")
    print("average classical communication cost is known to be exactly 1 bit.")
    print("\nFinal Equation for Communication (C):")
    print(f"C = {communication_cost} bit")

    # --- Resource 2: Non-signaling PR-Boxes ---
    chsh_classical_limit = 2.0
    chsh_quantum_limit = 2 * math.sqrt(2)
    chsh_pr_box_limit = 4.0

    # The amount of "non-locality" is the violation above the classical limit.
    quantum_non_locality = chsh_quantum_limit - chsh_classical_limit
    pr_box_non_locality = chsh_pr_box_limit - chsh_classical_limit

    # The cost is the ratio of the desired non-locality to the resource's non-locality.
    pr_box_cost = quantum_non_locality / pr_box_non_locality

    print("\n" + "="*60)
    print("Resource 2: Non-signaling PR-Boxes")
    print("="*60)
    print("This cost is calculated by comparing the violation of the CHSH inequality.")
    print(f"The classical LHV limit for the CHSH value is: {chsh_classical_limit}")
    print(f"The maximal quantum value (Tsirelson's bound) is: {chsh_quantum_limit:.4f}")
    print(f"The value for a PR-Box is: {chsh_pr_box_limit}")
    
    print("\nThe cost is the ratio of 'quantum non-locality' to 'PR-Box non-locality'.")
    print("\nFinal Equation for PR-Box cost (P):")
    print(f"P = (Quantum Limit - Classical Limit) / (PR-Box Limit - Classical Limit)")
    # Printing each number in the final equation
    print(f"P = ({chsh_quantum_limit:.4f} - {chsh_classical_limit}) / ({chsh_pr_box_limit} - {chsh_classical_limit})")
    print(f"P = {quantum_non_locality:.4f} / {pr_box_non_locality:.4f}")
    print(f"P = sqrt(2) - 1 ≈ {pr_box_cost:.4f} PR-Boxes (on average per trial)")

    # --- Resource 3: Shared Randomness ---
    # This is theoretically infinite for perfect simulation of all possible measurements.
    shared_randomness_cost = float('inf')
    print("\n" + "="*60)
    print("Resource 3: Shared Randomness (The 'Local Hidden Variable')")
    print("="*60)
    print("Any LHV model, by definition, requires a shared random variable.")
    print("To perfectly simulate the correlations for a continuous set of all possible")
    print("measurement settings, a source of infinite shared randomness is required.")
    
    print("\nFinal Equation for Shared Randomness (R):")
    print(f"R = {shared_randomness_cost}")


if __name__ == '__main__':
    calculate_simulation_resources()
