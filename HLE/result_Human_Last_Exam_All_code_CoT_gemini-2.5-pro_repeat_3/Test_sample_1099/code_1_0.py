import numpy as np
from scipy.integrate import quad

def simulate_singlet_cost():
    """
    Calculates the minimal average resources to simulate singlet state correlations.

    The problem is to find the average amount of non-signaling PR-Boxes and
    classical communication needed for a Local Hidden Variable (LHV) model to
    reproduce the correlations of a singlet state.

    1. The Model: We use the model by Brunner & Skrzypczyk. For a given angle 'theta'
       between Alice's and Bob's measurement settings, the quantum correlations can be
       simulated by a probabilistic mixture of a local model and a PR-box.

    2. PR-Box Cost: The probability of needing the PR-box for a given angle 'theta' is:
       p_PR(theta) = (1 - cos(theta)) / 2

    3. Communication Cost: This specific simulation model cleverly uses the PR-box
       to handle all the non-local aspects, and thus requires 0 bits of communication.

    4. Averaging: To find the average cost, we average over all possible relative
       measurement directions. The probability distribution for the angle 'theta' between
       two random vectors on a sphere is P(theta) = sin(theta) / 2.

    5. The Calculation: We need to compute the definite integral of the product of the
       PR-box probability and the angle distribution from 0 to pi.
       Average_PR_Cost = Integral from 0 to pi of [p_PR(theta) * P(theta)] d(theta)
    """

    # The integrand is p_PR(theta) * P(theta)
    # p_PR(theta) = (1 - cos(theta)) / 2
    # P(theta) = sin(theta) / 2
    integrand = lambda theta: ((1 - np.cos(theta)) / 2) * (np.sin(theta) / 2)

    # Perform the numerical integration from 0 to pi
    avg_pr_box_cost, _ = quad(integrand, 0, np.pi)

    # The communication cost in this model is zero
    avg_communication_cost = 0.0
    
    # --- Output ---
    print("To simulate the correlations of a singlet state with an LHV model, the minimal average resources required are:")
    print("-" * 80)
    
    # We are printing the components of the final calculation as requested
    print(f"The equation for the average PR-Box cost is:")
    print("Integral from 0 to pi of [ (1 - cos(theta))/2 * sin(theta)/2 ] d(theta)")
    
    print("\nFinal calculated average resources:")
    print(f"Average PR-Box units required = {avg_pr_box_cost:.4f}")
    print(f"Average bits of communication required = {avg_communication_cost:.4f}")
    print("-" * 80)
    
    # The final answer is the average number of PR-boxes required.
    # The user asked for "the minimal of ressources", which is a pair (PR-Box, Comm).
    # Since communication is 0, we can represent the key resource cost as a single number.
    
    # Return the answer in the required format.
    print(f"<<<{avg_pr_box_cost}>>>")


if __name__ == '__main__':
    simulate_singlet_cost()