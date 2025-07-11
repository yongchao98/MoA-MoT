import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a selfish mining attack.

    Args:
        beta (float): The adversary's fraction of the total mining power (0 < beta < 1).
        p (float): The probability that honest miners choose the adversary's block in a tie (0 <= p <= 1).
    """
    if not (0 < beta < 1):
        if beta == 0:
            print("With beta = 0, there is no adversary. Chain quality is 1.")
        elif beta >= 1:
            print("With beta >= 1, the adversary controls the chain. Chain quality is 0.")
        else:
            print("Error: beta must be between 0 and 1.", file=sys.stderr)
        return

    if not (0 <= p <= 1):
        print("Error: p must be between 0 and 1.", file=sys.stderr)
        return

    # Based on the Markov chain analysis, the final formula for chain quality (CQ) is:
    # CQ = ( (1-β)² * [1 + 2β(1-β)(1-p)] ) / ( 1 - β² + β³ )

    # Calculate the numerator parts
    # Part 1: (1 - beta)^2
    num_part1 = (1 - beta)**2
    # Part 2: 1 + 2 * beta * (1 - beta) * (1 - p)
    num_part2 = 1 + 2 * beta * (1 - beta) * (1 - p)
    numerator = num_part1 * num_part2

    # Calculate the denominator
    denominator = 1 - beta**2 + beta**3

    # Calculate the final chain quality
    chain_quality = numerator / denominator

    # Print the equation with all the numbers filled in
    print("Inputs:")
    print(f"  Adversary mining power (β) = {beta}")
    print(f"  Honest miner tie-breaking probability (p) = {p}\n")

    print("Formula for Chain Quality (CQ):")
    print("CQ = ( (1 - β)² * [1 + 2 * β * (1 - β) * (1 - p)] ) / (1 - β² + β³)\n")

    print("Calculation Steps:")
    print(f"Numerator = (1 - {beta})² * [1 + 2 * {beta} * (1 - {beta}) * (1 - {p})]")
    print(f"Numerator = {num_part1:.4f} * [1 + {2 * beta * (1 - beta) * (1 - p):.4f}]")
    print(f"Numerator = {num_part1:.4f} * {num_part2:.4f}")
    print(f"Numerator = {numerator:.4f}\n")

    print(f"Denominator = 1 - {beta}² + {beta}³")
    print(f"Denominator = 1 - {beta**2:.4f} + {beta**3:.4f}")
    print(f"Denominator = {denominator:.4f}\n")

    print("Final Result:")
    print(f"CQ = {numerator:.4f} / {denominator:.4f}")
    print(f"CQ = {chain_quality:.4f}")
    
    # Return the final numerical value for programmatic use
    return chain_quality


if __name__ == '__main__':
    # --- Parameters ---
    # β: Adversary's fraction of mining power.
    # p: Probability honest miners choose the adversary's block in a tie.
    beta_power = 0.35
    tie_prob = 0.25
    
    result = calculate_chain_quality(beta_power, tie_prob)
    # The final answer format is requested at the very end.
    # We will format the result to be included there.
    # Note: The result printed above is for user's understanding, the final tagged answer will be separate.
    # <<< f-string formatting will be resolved at execution time for the final value.
    print(f"\n<<<{result:.16f}>>>")
