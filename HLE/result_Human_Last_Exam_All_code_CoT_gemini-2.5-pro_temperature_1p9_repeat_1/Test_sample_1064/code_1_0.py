import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality for a selfish mining scenario.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability that honest miners choose the adversary's block in a tie (0 <= p <= 1).
    
    Returns:
        float: The expected chain quality.
    """
    if not (0 < beta < 1):
        print("Error: beta must be between 0 and 1.", file=sys.stderr)
        return None
    if not (0 <= p <= 1):
        print("Error: p must be between 0 and 1.", file=sys.stderr)
        return None

    alpha = 1.0 - beta

    # Based on the derivation, the expected chain quality is E_H / (E_H + E_A),
    # which simplifies to the following expression:
    
    # Numerator of the chain quality formula
    # (1-beta)^2 * (1 + beta*(1-beta)*(2-p))
    numerator = (alpha**2) * (1 + beta * alpha * (2 - p))
    
    # Denominator of the chain quality formula
    # 1 - beta^2 + beta^3
    denominator = 1 - beta**2 + beta**3

    if denominator == 0:
        print("Error: Denominator is zero, cannot calculate chain quality.", file=sys.stderr)
        return None
        
    chain_quality = numerator / denominator
    return numerator, denominator, chain_quality

def main():
    """
    Main function to execute the calculation and print the results.
    """
    # --- You can change these values ---
    # beta: adversary's mining power portion (e.g., 0.33 for 33%)
    beta = 0.33
    # p: probability of honest miners choosing the adversary's block in a tie
    p = 0.5
    # ------------------------------------

    print(f"Calculating Chain Quality for beta = {beta} and p = {p}\n")

    print("The final formula for Chain Quality (CQ) is:")
    print("CQ = [ (1-β)² * (1 + β*(1-β)*(2-p)) ] / [ 1 - β² + β³ ]\n")

    result = calculate_chain_quality(beta, p)
    
    if result is not None:
        numerator, denominator, quality = result
        print("Substituting the values into the formula:")
        print(f"Numerator   = (1-{beta})² * (1 + {beta}*(1-{beta})*(2-{p})) = {numerator}")
        print(f"Denominator = 1 - {beta}² + {beta}³ = {denominator}")
        print(f"\nExpected Chain Quality = {numerator} / {denominator}")
        print(f"                       = {quality}")
        print("\nThis means that in the long run, approximately {:.2f}% of the blocks in the longest chain will be honest blocks.".format(quality * 100))
        
        # The final answer in the requested format
        print(f"\n<<<{quality}>>>")

if __name__ == "__main__":
    main()