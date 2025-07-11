import math

def calculate_simulation_resources():
    """
    Calculates the minimal resources (communication or PR-Boxes) needed to
    simulate the correlations of measurements on a singlet state with an LHV model.
    """

    # --- Resource 1: Classical Communication ---
    # According to the Toner-Bacon protocol (2003), the average classical
    # communication cost to perfectly simulate the correlations from any
    # projective measurement on a singlet state is exactly 1 bit.
    communication_cost = 1.0

    # --- Resource 2: Non-signaling PR-Box ---
    # The cost in terms of PR-Boxes is calculated by finding the fraction of a
    # PR-Box needed to achieve the maximum quantum violation of the CHSH inequality
    # (the Tsirelson bound).

    # CHSH inequality bounds for different models:
    classical_limit = 2.0
    quantum_limit = 2 * math.sqrt(2)
    pr_box_limit = 4.0

    # The cost is the fraction 'p' in a mixture p*PR + (1-p)*LocalBox
    # that yields the quantum limit.
    # Formula: p = (Quantum - Classical) / (PR-Box - Classical)
    pr_box_cost = (quantum_limit - classical_limit) / (pr_box_limit - classical_limit)

    # --- Print the results ---
    print("To supplement a Local Hidden Variable (LHV) model to simulate the correlations of a singlet state, the average minimal resource cost is:")

    print("\n--- OPTION 1: Using Classical Communication ---")
    print(f"The model requires an average of {communication_cost} bit of classical communication.")
    print("Final Equation:")
    print(f"C_comm = {communication_cost}")

    print("\n--- OPTION 2: Using a Non-Signaling PR-Box ---")
    print("The model requires a fraction of a PR-Box resource per trial.")
    print("Final Equation:")
    print(f"C_prbox = (Quantum_Limit - Classical_Limit) / (PR_Box_Limit - Classical_Limit)")
    # We output the numbers in the final equation as requested
    print(f"C_prbox = (2*sqrt(2) - {classical_limit}) / ({pr_box_limit} - {classical_limit})")
    print(f"Result: C_prbox = ({quantum_limit:.4f} - {classical_limit}) / ({pr_box_limit} - {classical_limit}) = {pr_box_cost:.4f}")

if __name__ == '__main__':
    calculate_simulation_resources()
