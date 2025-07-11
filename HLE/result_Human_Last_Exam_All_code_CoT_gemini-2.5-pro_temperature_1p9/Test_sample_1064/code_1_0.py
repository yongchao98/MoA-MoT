import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a selfish mining attack.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability an honest miner chooses the adversary's block in a tie (0 <= p <= 1).
        
    Returns:
        float: The expected chain quality.
    """
    # Input validation
    if not (0 < beta < 1):
        print("Error: beta must be between 0 and 1.", file=sys.stderr)
        return None
    if not (0 <= p <= 1):
        print("Error: p must be between 0 and 1.", file=sys.stderr)
        return None

    alpha = 1.0 - beta

    # The final formula for chain quality Q is derived from a Markov chain model:
    # Q = (alpha^2 * (1 + alpha*beta*(2-p))) / (1 - beta^2 + beta^3)

    # Calculate the numerator of the formula
    # Numerator = (1-β)^2 * (1 + (1-β)β(2-p))
    num_term1 = alpha**2
    num_term2 = 1 + alpha * beta * (2 - p)
    numerator = num_term1 * num_term2
    
    # Calculate the denominator of the formula
    # Denominator = 1 - β^2 + β^3
    denominator = 1 - beta**2 + beta**3

    if denominator == 0:
        print("Error: Denominator is zero. Cannot calculate chain quality.", file=sys.stderr)
        return None

    # Calculate the final chain quality
    chain_quality = numerator / denominator
    
    # Print out the components of the final equation as requested
    print("--- Calculation Details ---")
    print(f"Inputs: beta = {beta}, p = {p}")
    print(f"alpha (1-beta) = {alpha:.4f}\n")
    
    print("Numerator Formula: (alpha^2) * (1 + alpha * beta * (2 - p))")
    print(f"  alpha^2 = {num_term1:.4f}")
    print(f"  1 + alpha*beta*(2-p) = {num_term2:.4f}")
    print(f"Final Numerator Value = {numerator:.4f}\n")

    print("Denominator Formula: 1 - beta^2 + beta^3")
    print(f"  1 - beta^2 = {1 - beta**2:.4f}")
    print(f"  beta^3 = {beta**3:.4f}")
    print(f"Final Denominator Value = {denominator:.4f}\n")
    
    print("--- Final Result ---")
    print(f"Expected Chain Quality = Numerator / Denominator")
    print(f"                       = {numerator:.4f} / {denominator:.4f}")
    print(f"                       = {chain_quality:.6f}")

    return chain_quality

if __name__ == '__main__':
    # Example usage:
    # Set the adversary's mining power (beta)
    # For selfish mining to be profitable, beta usually needs to be > 1/4 or 1/3 depending on network conditions
    adversary_power_beta = 0.35
    
    # Set the probability p for tie-breaking
    # p=0 means honest miners never adopt the selfish miner's block in a tie
    # p=1 means honest miners always adopt the selfish miner's block in a tie
    # p=0.5 means random choice
    tie_breaking_prob_p = 0.5
    
    calculate_chain_quality(adversary_power_beta, tie_breaking_prob_p)