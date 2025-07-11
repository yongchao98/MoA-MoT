import math

def calculate_simulation_resources():
    """
    Calculates and prints the minimal average resources needed to simulate
    the correlations of a singlet quantum state with an LHV model.
    """

    # --- Resource 1: Classical Communication (C) ---
    # Based on the Toner-Bacon protocol (2003), a perfect simulation can be
    # achieved by supplementing an LHV model with a specific amount of
    # one-way classical communication.
    communication_cost_bits = 1.0

    # --- Resource 2: Equivalent PR-Box non-locality (N_PR) ---
    # This measures the non-local strength of a singlet state relative to a
    # Popescu-Rohrlich (PR) box by comparing CHSH inequality violations.

    # CHSH value for a Local Hidden Variable (LHV) model (classical limit)
    chsh_classical = 2.0

    # CHSH value for a singlet quantum state (Tsirelson's bound)
    chsh_quantum = 2 * math.sqrt(2)

    # CHSH value for a non-signaling PR-Box (algebraic maximum)
    chsh_pr_box = 4.0

    # The "non-local advantage" is the amount the score exceeds the classical limit.
    advantage_quantum = chsh_quantum - chsh_classical
    advantage_pr_box = chsh_pr_box - chsh_classical

    # The required fraction of PR-Box non-locality is the ratio of the advantages.
    n_pr = advantage_quantum / advantage_pr_box

    print("The minimal average resources required to simulate singlet state correlations with an LHV model are:")
    print("-" * 80)
    print(f"1. Classical Communication Cost (C): {communication_cost_bits} bit")
    print("   Explanation: An exact simulation is possible using 1 bit of communication (Toner-Bacon protocol).")
    print("-" * 80)
    print("2. Equivalent PR-Box Cost (N_PR):")
    print("   Explanation: This quantifies the singlet's non-local strength relative to a PR-Box.")
    print("\n   The calculation uses the maximal CHSH inequality scores:")
    print(f"   - Classical Limit: {chsh_classical}")
    print(f"   - Quantum Limit (2 * sqrt(2)): {chsh_quantum:.8f}")
    print(f"   - PR-Box Limit: {chsh_pr_box}")

    print("\n   The final calculation for the PR-Box fraction is:")
    # As requested, printing each number in the final equation
    print(f"   N_PR = ( {chsh_quantum:.8f} - {chsh_classical} ) / ( {chsh_pr_box} - {chsh_classical} )")
    print(f"   N_PR = {n_pr:.8f}")
    print("\n   This indicates that a singlet state carries approx. 41.4% of the non-local power of a PR-Box.")

if __name__ == '__main__':
    calculate_simulation_resources()
