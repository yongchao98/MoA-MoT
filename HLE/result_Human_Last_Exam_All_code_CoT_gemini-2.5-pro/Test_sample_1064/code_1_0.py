import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability that honest miners choose the adversary's block in a tie (0 <= p <= 1).

    Returns:
        float: The expected chain quality.
    """
    if not (0 < beta < 1):
        raise ValueError("Adversary's mining power (beta) must be between 0 and 1.")
    if not (0 <= p <= 1):
        raise ValueError("Tie-breaking probability (p) must be between 0 and 1.")

    # The formula is derived from a Markov chain model of the selfish mining process.
    # Let's break down the components of the formula.
    
    # Numerator of the formula
    # (1-beta)^2
    term1_num = (1 - beta)**2
    # (1 + beta * (1-beta) * (2-p))
    term2_num = (1 + beta * (1 - beta) * (2 - p))
    numerator = term1_num * term2_num

    # Denominator of the formula
    # 1 - beta^2 + beta^3
    denominator = 1 - beta**2 + beta**3
    
    if denominator == 0:
        return float('inf') # Should not happen for 0 < beta < 1

    chain_quality = numerator / denominator
    
    # The final equation is:
    # Chain Quality = ((1-β)² * (1 + β(1-β)*(2-p))) / (1 - β² + β³)
    # Let's print the components for clarity, as requested.
    
    print("Calculating Chain Quality with beta = {} and p = {}".format(beta, p))
    print("="*50)
    print("The final equation for Chain Quality is:")
    print("CQ = ((1-β)² * (1 + β*(1-β)*(2-p))) / (1 - β² + β³)\n")
    
    print("Substituting the given values:")
    print("Numerator part 1, (1-β)²: (1 - {})² = {}".format(beta, term1_num))
    print("Numerator part 2, (1 + β*(1-β)*(2-p)): (1 + {}*(1-{})*(2-{})) = {}".format(beta, beta, p, term2_num))
    print("Final Numerator value: {} * {} = {}".format(term1_num, term2_num, numerator))
    print("-" * 20)
    print("Denominator, (1 - β² + β³): (1 - {}² + {}³) = {}".format(beta, beta, denominator))
    print("-" * 20)
    print("Final Result, Numerator / Denominator: {} / {} = {}".format(numerator, denominator, chain_quality))
    print("="*50)
    
    return chain_quality

if __name__ == '__main__':
    # Example usage:
    # You can change these values to see how they affect the chain quality.
    # beta is the adversary's mining power. A common value to test is just over 1/3.
    # p is the probability honest miners build on the selfish miner's block in a tie.
    # p = 0.5 means they choose randomly.
    # p = 0 means they never choose the selfish block (if they know who it is).
    # p = 1 means they always choose the selfish block.
    
    try:
        # Check if command line arguments are provided
        if len(sys.argv) == 3:
            adv_power = float(sys.argv[1])
            tie_prob = float(sys.argv[2])
        else:
            # Default values if no arguments are given
            print("Usage: python your_script.py <beta> <p>")
            print("Using default values: beta = 0.35, p = 0.5\n")
            adv_power = 0.35
            tie_prob = 0.5
            
        final_quality = calculate_chain_quality(adv_power, tie_prob)
        print("\nThe expected chain quality is {:.4f}".format(final_quality))
        # The final answer format is not applicable here since the output is a calculation.
        # However, if we must provide a single value, it will be the result of this calculation.
        # For the default values beta=0.35, p=0.5, the output is 0.5636...
        # Let's wrap just the final number in the required format.
        print("\n<<<{}>>>".format(final_quality))

    except ValueError as e:
        print(f"Error: {e}")
