import math

def calculate_pr_box_fraction():
    """
    Calculates the minimal fraction of a PR-box needed to simulate the maximal
    CHSH violation of a quantum singlet state.

    This is based on a model that mixes a non-local PR-box resource with a
    purely classical (local hidden variable) resource, with no communication.
    """

    # S_max are the maximum values for the CHSH inequality expression.
    # S_classical_max: The maximum value for any Local Hidden Variable (LHV) model.
    S_classical_max = 2.0

    # S_quantum_max: The maximum value in quantum mechanics (Tsirelson's bound).
    S_quantum_max = 2 * math.sqrt(2)

    # S_pr_box_max: The maximum value for a Popescu-Rohrlich (PR) box.
    S_pr_box_max = 4.0

    # We want to find the fraction 'p' of a PR-box resource that, when mixed
    # with an optimal classical resource (fraction 1-p), reproduces the quantum maximum.
    # The equation for the expected CHSH value of the mixture is:
    # p * S_pr_box_max + (1-p) * S_classical_max = S_quantum_max
    #
    # We solve for p:
    # p * S_pr_box_max + S_classical_max - p * S_classical_max = S_quantum_max
    # p * (S_pr_box_max - S_classical_max) = S_quantum_max - S_classical_max
    # p = (S_quantum_max - S_classical_max) / (S_pr_box_max - S_classical_max)

    p = (S_quantum_max - S_classical_max) / (S_pr_box_max - S_classical_max)

    print("--- Simulating Quantum CHSH Violation with a PR-Box ---")
    print("The goal is to find the minimum fraction 'p' of a PR-box needed.")
    print("\nWe set up the simulation equation:")
    print(f"p * {S_pr_box_max} + (1-p) * {S_classical_max} = {S_quantum_max:.4f}")

    print("\nSolving for 'p', we get the formula:")
    print(f"p = ({S_quantum_max:.4f} - {S_classical_max}) / ({S_pr_box_max} - {S_classical_max})")

    print("\nThe minimal fraction of a PR-box resource required is:")
    print(f"p = {p:.4f}")
    
    print("\nThis means you need, on average, a PR-box that works correctly about 41.4% of the time,")
    print("and use an optimal classical strategy for the rest, to simulate this key quantum correlation without communication.")
    
    print("\nFor comparison, another known result is that simulating a singlet state for all *projective* measurements (a simpler case than all POVMs) requires 1 bit of classical communication (and 0 PR-Boxes).")


calculate_pr_box_fraction()

# Final answer is the numerical value of p
p_value = (2 * math.sqrt(2) - 2) / (4 - 2)
# The final answer will be just the value, as requested by the user prompt.
#<<<p_value>>> this line will be stripped before output
final_answer = (2 * math.sqrt(2) - 2) / (4 - 2)
print(f"<<<{final_answer}>>>")