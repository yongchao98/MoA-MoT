import sys

def solve_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The adversary's portion of mining power (0 < beta < 1).
        p (float): The probability of honest miners choosing the adversary's block in a tie (0 <= p <= 1).
    """
    if not (0 < beta < 1):
        print("Error: Adversary's mining power (β) must be between 0 and 1.", file=sys.stderr)
        return
    if not (0 <= p <= 1):
        print("Error: Tie-breaking probability (p) must be between 0 and 1.", file=sys.stderr)
        return

    alpha = 1.0 - beta

    # Calculate numerator and denominator of the chain quality formula
    # Formula: Q = (α^2 + α^3*β*(2-p)) / (α^2 + α^2*β + β)
    numerator = alpha**2 + alpha**3 * beta * (2 - p)
    denominator = alpha**2 + alpha**2 * beta + beta

    # Calculate the final chain quality
    chain_quality = numerator / denominator

    # Print the explanation and results
    print("This script calculates the expected chain quality in a Bitcoin system with a selfish miner.")
    print(f"Given parameters:\n- Adversary's mining power (β): {beta:.4f}\n- Honest miner tie-breaking probability (p): {p:.4f}\n")
    print(f"Derived parameter:\n- Honest mining power (α = 1 - β): {alpha:.4f}\n")

    print("The formula for the expected chain quality (Q) is:")
    print("Q = (α^2 + α^3*β*(2-p)) / (α^2 + α^2*β + β)\n")

    print("Plugging in the values for the equation:")
    print(f"Numerator   = ({alpha:.4f})^2 + ({alpha:.4f})^3 * {beta:.4f} * (2 - {p:.4f})")
    numerator_part1 = alpha**2
    numerator_part2 = alpha**3 * beta * (2 - p)
    print(f"            = {numerator_part1:.4f} + {numerator_part2:.4f}")
    print(f"            = {numerator:.4f}\n")

    print(f"Denominator = ({alpha:.4f})^2 + ({alpha:.4f})^2 * {beta:.4f} + {beta:.4f}")
    denominator_part1 = alpha**2
    denominator_part2 = alpha**2 * beta
    denominator_part3 = beta
    print(f"            = {denominator_part1:.4f} + {denominator_part2:.4f} + {denominator_part3:.4f}")
    print(f"            = {denominator:.4f}\n")

    print("--- Final Result ---")
    print(f"Expected Chain Quality (Q) = {numerator:.4f} / {denominator:.4f} = {chain_quality:.4f}")

    # Return the final numeric answer for the specified format.
    # Note: The '<<<' format is for automated processing and should contain only the final numerical value.
    print(f"\n<<<{chain_quality:.4f}>>>")


if __name__ == '__main__':
    # --- Parameters ---
    # β: The portion of mining power controlled by the adversary.
    # A common value used for analysis is 1/3.
    adversary_power_beta = 1.0/3.0

    # p: The probability that honest miners choose the adversary's block in a tie.
    # p = 0.5 means the network is split evenly.
    tie_breaking_prob_p = 0.5

    solve_chain_quality(adversary_power_beta, tie_breaking_prob_p)