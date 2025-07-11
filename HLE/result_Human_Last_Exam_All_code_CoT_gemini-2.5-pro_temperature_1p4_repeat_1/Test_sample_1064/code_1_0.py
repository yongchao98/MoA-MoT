import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a specific selfish mining model.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability an honest miner chooses the adversary's block in a tie (0 <= p <= 1).

    Returns:
        float: The expected chain quality. Returns -1.0 if inputs are invalid.
    """
    if not (0 < beta < 1):
        print("Error: beta must be between 0 and 1.", file=sys.stderr)
        return -1.0
    if not (0 <= p <= 1):
        print("Error: p must be between 0 and 1.", file=sys.stderr)
        return -1.0
        
    # The derived formula is valid for beta < 0.5. For beta >= 0.5,
    # the adversary's lead tends to grow infinitely, making the chain
    # quality tend towards 0.
    if beta >= 0.5:
        print("For beta >= 0.5, the chain quality approaches 0 as the adversary's chain will grow indefinitely.")
        numerator = 0
        denominator = 1
        quality = 0
    else:
        # Calculate the expected number of honest blocks in a cycle.
        # h_0 = 1 - beta + beta * (1-beta)^2 * (2-p)
        term1_h = 1 - beta
        term2_h = beta * ((1 - beta) ** 2) * (2 - p)
        numerator = term1_h + term2_h

        # Calculate the expected total number of blocks in a cycle.
        # L_0 = (1 - beta^2 + beta^3) / (1-beta)
        denominator_num = 1 - (beta ** 2) + (beta ** 3)
        denominator_den = 1 - beta
        denominator = denominator_num / denominator_den
        
        # Chain Quality is the ratio of the two.
        quality = numerator / denominator

    print("For beta = {} and p = {}:".format(beta, p))
    print("Final Equation: Chain Quality = Numerator / Denominator")
    print("-------------------------------------------------------")
    print(f"Numerator (Expected Honest Blocks per Cycle): {numerator}")
    print(f"Denominator (Expected Total Blocks per Cycle): {denominator}")
    print("-------------------------------------------------------")
    print(f"Expected Chain Quality: {quality}")

if __name__ == '__main__':
    # You can change these values to explore different scenarios.
    # beta is the adversary's mining power.
    # p is the probability honest miners choose the adversary's block in a tie.
    
    beta_power = 0.3
    tie_prob = 0.5

    calculate_chain_quality(beta_power, tie_prob)

    print("\nExample for the original selfish mining case (p=0):")
    calculate_chain_quality(0.3, 0)