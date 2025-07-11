import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 0.5).
        p (float): The probability that an honest miner chooses the adversary's block
                   in a tie (0 <= p <= 1).

    Returns:
        float: The expected chain quality.
    """
    if not 0 < beta < 0.5:
        print("Error: The adversary's mining power 'beta' must be between 0 and 0.5 for the chain to be stable.", file=sys.stderr)
        # For beta >= 0.5, the adversary's lead grows infinitely, chain quality approaches 0.
        return 0.0

    if not 0 <= p <= 1:
        print("Error: The probability 'p' must be between 0 and 1.", file=sys.stderr)
        return -1.0

    print(f"Calculating chain quality for beta = {beta} and p = {p}\n")
    
    b2 = beta ** 2
    b3 = beta ** 3

    # Numerator terms of the formula: 1 - p*β + (3p-7)β^2 + (6-2p)β^3
    num_term1 = -p * beta
    num_term2 = (3 * p - 7) * b2
    num_term3 = (6 - 2 * p) * b3
    numerator = 1 + num_term1 + num_term2 + num_term3

    # Denominator terms of the formula: 1 - 3β^2 - β^3
    den_term1 = -3 * b2
    den_term2 = -b3
    denominator = 1 + den_term1 + den_term2
    
    if denominator == 0:
        print("Error: Calculation resulted in division by zero.", file=sys.stderr)
        return -1.0
        
    chain_quality = numerator / denominator

    # Output the details of the calculation as per the final equation
    print("The final equation for chain quality (CQ) is:")
    print("CQ = (1 - p*β + (3p-7)β^2 + (6-2p)β^3) / (1 - 3β^2 - β^3)\n")
    
    print("Substituting the given values:")
    print("Numerator calculation:")
    print(f"  1 + ({num_term1:.4f}) + ({num_term2:.4f}) + ({num_term3:.4f}) = {numerator:.4f}")

    print("\nDenominator calculation:")
    print(f"  1 + ({den_term1:.4f}) + ({den_term2:.4f}) = {denominator:.4f}")

    print("\nFinal Calculation:")
    print(f"  CQ = {numerator:.4f} / {denominator:.4f}")
    
    return chain_quality

if __name__ == '__main__':
    # --- Parameters ---
    # β: Adversary's mining power. Must be < 0.5 for a stable result.
    adversary_power_beta = 0.25
    # p: Probability honest miners choose the adversary's block in a tie.
    tie_breaking_prob_p = 0.5

    # --- Calculation ---
    quality = calculate_chain_quality(adversary_power_beta, tie_breaking_prob_p)
    
    if quality >= 0:
        print(f"\nExpected Chain Quality: {quality:.4f}")
        print(f"<<<{quality:.4f}>>>")
