import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The adversary's portion of the mining power (0 <= beta < 1).
        p (float): The probability that an honest miner chooses the adversary's 
                   block in a tie (0 <= p <= 1).
    """
    # --- Input Validation ---
    if not (0 <= beta < 1):
        print("Error: beta (adversary's mining power) must be between 0 and 1.", file=sys.stderr)
        return
    if not (0 <= p <= 1):
        print("Error: p (tie-breaking probability) must be between 0 and 1.", file=sys.stderr)
        return

    # --- Calculation ---
    # If the adversary's power is 50% or more, their lead grows indefinitely.
    # The main chain will consist only of their blocks in the long run.
    if beta >= 0.5:
        chain_quality = 0.0
        print(f"Given beta = {beta}, which is >= 0.5, the selfish miner's lead will grow indefinitely.")
        print("The long-term chain will consist exclusively of the adversary's blocks.")
        print("\nExpected Chain Quality = 0.0")
        print(f"<<<{chain_quality}>>>")
        return

    # For beta < 0.5, we use the derived formula for chain quality Q.
    # Q = ((1-beta)*(1 + 2*beta - beta*p) * (1-2*beta)) / (1 - beta - beta^2)
    
    # Calculate each part of the formula
    num_part1 = 1 - beta
    num_part2 = 1 + 2 * beta - beta * p
    num_part3 = 1 - 2 * beta
    numerator = num_part1 * num_part2 * num_part3
    
    den_part1 = 1
    den_part2 = -beta
    den_part3 = -beta**2
    denominator = den_part1 + den_part2 + den_part3
    
    if denominator == 0:
        print("Error: Denominator is zero. Cannot calculate chain quality.", file=sys.stderr)
        return
        
    chain_quality = numerator / denominator
    
    # --- Output ---
    print("The formula for Expected Chain Quality (Q) for beta < 0.5 is:")
    print("Q = [(1 - beta) * (1 + 2*beta - beta*p) * (1 - 2*beta)] / [1 - beta - beta^2]")
    print("\nPlugging in the values beta = {} and p = {}:".format(beta, p))
    
    print("\nNumerator = (1 - {}) * (1 + 2*{} - {}*{}) * (1 - 2*{})".format(beta, beta, beta, p, beta))
    print("          = ({}) * ({}) * ({})".format(num_part1, num_part2, num_part3))
    print("          = {}".format(numerator))
    
    print("\nDenominator = 1 - {} - {}^2".format(beta, beta))
    print("            = {} + ({}) + ({})".format(den_part1, den_part2, den_part3))
    print("            = {}".format(denominator))
    
    print("\nQ = Numerator / Denominator")
    print("Q = {} / {}".format(numerator, denominator))
    print("\nFinal Expected Chain Quality = {}".format(chain_quality))
    print(f"<<<{chain_quality}>>>")


if __name__ == '__main__':
    # --- Parameters ---
    # You can change these values to see the effect on chain quality.
    
    # beta: The adversary's portion of the mining power.
    # Example: 0.3 means 30% of the total hash power.
    adversary_power_beta = 0.3

    # p: The probability that an honest miner chooses the adversary's block in a tie.
    # p = 0.5 means honest miners are indifferent.
    # p = 0 means honest miners always prefer their own block.
    # p = 1 means honest miners are forced to build on the adversary's block.
    tie_breaking_prob_p = 0.5
    
    calculate_chain_quality(adversary_power_beta, tie_breaking_prob_p)