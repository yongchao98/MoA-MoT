import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability that an honest miner chooses the adversary's block
                   in a tie (0 <= p <= 1).

    Returns:
        float: The expected chain quality.
    """
    if not (0 < beta < 1):
        raise ValueError("Adversary's mining power (beta) must be between 0 and 1.")
    if not (0 <= p <= 1):
        raise ValueError("The tie-breaking probability (p) must be between 0 and 1.")

    alpha = 1.0 - beta

    # Based on the state transition analysis, the rate of blocks being added
    # to the chain is alpha. The rate of honest blocks being added is R_H.
    # The chain quality Q = R_H / alpha.

    # We first calculate pi_0, the steady-state probability of being in state 0 (no lead).
    # pi_0 = alpha / (alpha + alpha^2 * beta + beta)
    pi_0_denominator = alpha + (alpha**2 * beta) + beta
    pi_0 = alpha / pi_0_denominator

    # The chain quality Q is given by the formula: Q = pi_0 * (1 + alpha * beta * (2 - p))
    # We calculate the 'factor' part of the equation first.
    quality_factor = 1 + alpha * beta * (2 - p)

    # The final chain quality is the product of pi_0 and the factor.
    chain_quality = pi_0 * quality_factor

    print("### Calculation Steps ###")
    print(f"Given Parameters:")
    print(f"  Adversary mining power (β): {beta}")
    print(f"  Honest miner tie-choice probability (p): {p}\n")

    print(f"Intermediate Values:")
    print(f"  Honest mining power (α = 1 - β): {alpha:.4f}")
    
    # Printing the components of the final formula as requested.
    # The formula is: chain_quality = pi_0 * (1 + alpha * beta * (2-p))
    print(f"  pi_0 (probability of state 0): {pi_0:.4f}")
    print(f"  (1 + alpha * beta * (2-p)): {quality_factor:.4f}")

    print("\n---")
    print("Final Result:")
    print(f"The expected chain quality is: {chain_quality:.6f}")
    
    return chain_quality

if __name__ == '__main__':
    # Example Usage:
    # You can change these values to see the effect on chain quality.
    # beta is the adversary's portion of mining power.
    # p is the probability an honest miner chooses the adversary's block in a tie.
    try:
        # Check if command line arguments are provided
        if len(sys.argv) == 3:
            adv_power_beta = float(sys.argv[1])
            tie_break_p = float(sys.argv[2])
            print(f"Running with provided arguments: β = {adv_power_beta}, p = {tie_break_p}\n")
        else:
            adv_power_beta = 0.333
            tie_break_p = 0.5
            print(f"Running with default example values: β = {adv_power_beta}, p = {tie_break_p}")
            print("You can also provide values from the command line, e.g.,")
            print(f"python {sys.argv[0]} 0.4 0.25\n")
        
        final_quality = calculate_chain_quality(beta=adv_power_beta, p=tie_break_p)
        # Final answer formatted for parsing
        print(f"\n<<< {final_quality} >>>")

    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
