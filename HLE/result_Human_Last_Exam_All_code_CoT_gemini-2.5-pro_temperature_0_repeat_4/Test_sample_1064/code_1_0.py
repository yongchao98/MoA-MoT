import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The adversary's portion of the total mining power (0 < beta < 1).
        p (float): The probability that honest miners choose the adversary's block in a tie (0 <= p <= 1).
    """
    # Validate inputs
    if not (0 < beta < 1):
        print("Error: beta (adversary's mining power) must be between 0 and 1.", file=sys.stderr)
        return
    if not (0 <= p <= 1):
        print("Error: p (tie-breaking probability) must be between 0 and 1.", file=sys.stderr)
        return

    # The formula for chain quality (Q) is derived from a Markov chain model of selfish mining:
    # Q = ( (1-β)² * (1 + β*(1-β)*(2-p)) ) / ( 1 - β² + β³ )

    # Calculate numerator and denominator for clarity
    numerator = ((1 - beta)**2) * (1 + beta * (1 - beta) * (2 - p))
    denominator = 1 - beta**2 + beta**3

    if denominator == 0:
        print("Error: Calculation resulted in division by zero.", file=sys.stderr)
        return

    # Calculate the final chain quality
    chain_quality = numerator / denominator

    # Print the results, showing the formula with the numbers substituted
    print(f"Given:")
    print(f"  Adversary's mining power (β) = {beta}")
    print(f"  Tie-breaking probability (p) = {p}\n")

    print(f"The formula for expected chain quality (Q) is:")
    print(f"  Q = ( (1-β)² * (1 + β*(1-β)*(2-p)) ) / ( 1 - β² + β³ )\n")

    print(f"Substituting the given values:")
    # Show each number in the final equation
    print(f"  Q = ( (1-{beta})² * (1 + {beta}*(1-{beta})*(2-{p})) ) / ( 1 - {beta}² + {beta}³ )")
    
    # Show intermediate calculations
    term1 = (1 - beta)**2
    term2 = 1 + beta * (1 - beta) * (2 - p)
    print(f"  Q = ( {term1:.4f} * {term2:.4f} ) / ( {denominator:.4f} )")
    print(f"  Q = {numerator:.4f} / {denominator:.4f}\n")

    print(f"Result:")
    print(f"  Expected Chain Quality = {chain_quality:.4f}")
    print(f"  This means we expect about {chain_quality:.2%} of the blocks in the longest chain to be from honest miners.")


if __name__ == '__main__':
    # Example values from academic literature on selfish mining
    # An adversary with 1/3 of the mining power
    adversary_power_beta = 1/3
    # A scenario where in case of a tie, the network is split 50/50
    tie_breaking_p = 0.5

    calculate_chain_quality(adversary_power_beta, tie_breaking_p)
