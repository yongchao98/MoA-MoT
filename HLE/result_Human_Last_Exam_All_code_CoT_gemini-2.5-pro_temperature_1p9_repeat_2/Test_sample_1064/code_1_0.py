import argparse

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a selfish mining attack.

    Args:
        beta (float): The fraction of mining power controlled by the adversary. 
                      Must be in the range [0, 1).
        p (float): The probability that an honest miner chooses the adversary's block 
                   in a tie. Must be in the range [0, 1].

    Returns:
        float: The expected chain quality.
    """
    if not 0 <= beta < 1:
        raise ValueError("Beta must be in the range [0, 1).")
    if not 0 <= p <= 1:
        raise ValueError("p must be in the range [0, 1].")

    # Let alpha be the honest mining power
    alpha = 1.0 - beta

    # Expected number of honest blocks added to the chain per "block event".
    # This is derived by summing the honest blocks from all possible outcomes
    # (honest discovery, or a selfish mining episode), weighted by their probabilities.
    # E_h = (1-β)*1 (honest block) + β*(...)(selfish mining outcomes)
    expected_honest_blocks = alpha + beta * (alpha**2) * (2 - p)

    # Expected total number of blocks added to the chain per "block event".
    # This is also a weighted sum over all outcomes.
    # The term for the selfish part (β*...) accounts for ties and the adversary
    # gaining a lead of 2 or more.
    if alpha == 0: # Technically beta cannot be 1, but for safety.
        # If beta is close to 1, this part dominates. As beta -> 1, chain quality -> 0.
        expected_total_blocks = float('inf')
    else:
        # Denominator from the full formula: (1 - beta^2 + beta^3) / (1 - beta)
        numerator = 1 - beta**2 + beta**3
        expected_total_blocks = numerator / alpha
    
    if expected_total_blocks == 0:
        return 0.0

    chain_quality = expected_honest_blocks / expected_total_blocks
    
    print("This script calculates the expected chain quality given the adversary's mining power (beta) and tie-breaking probability (p).")
    print("-" * 30)
    print(f"Inputs: beta = {beta}, p = {p}")
    print("-" * 30)
    print(f"The calculation is based on the formula: Chain Quality = E_h / E_l")
    print(f"Where E_h is the expected number of honest blocks per round, and E_l is the expected total blocks per round.")
    print(f"\nCalculated value for E_h: {expected_honest_blocks}")
    print(f"Calculated value for E_l: {expected_total_blocks}")
    print(f"\nFinal Equation:")
    print(f"{expected_honest_blocks} / {expected_total_blocks} = {chain_quality}")
    
    return chain_quality

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate expected chain quality under selfish mining.")
    parser.add_argument('--beta', type=float, default=0.3, help="Adversary's fraction of mining power (e.g., 0.3 for 30%%). Must be between 0 and 1.")
    parser.add_argument('--p', type=float, default=0.5, help="Probability honest miners choose the adversary's block in a tie (e.g., 0.5). Must be between 0 and 1.")
    args = parser.parse_args()
    
    calculate_chain_quality(args.beta, args.p)