import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a specific selfish mining model.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability of the network choosing the adversary's block
                   in a tie (0 <= p <= 1).
    """
    if not (0 < beta < 1):
        print("Error: Adversary's mining power (beta) must be between 0 and 1.", file=sys.stderr)
        return

    if not (0 <= p <= 1):
        print("Error: Tie-breaking probability (p) must be between 0 and 1.", file=sys.stderr)
        return

    alpha = 1 - beta

    # Numerator of the chain quality formula: α(1 + 2αβ - αβp)
    # This represents the weighted rate of honest blocks being added to the chain.
    numerator = alpha + 2 * (alpha**2) * beta - (alpha**2) * beta * p

    # Denominator of the chain quality formula: αβ + α + β/α
    # This represents the weighted rate of total (honest + adversary) blocks being added.
    denominator = (alpha * beta) + alpha + (beta / alpha)

    # Expected Chain Quality
    chain_quality = numerator / denominator

    # Output the equation with the numbers used in the calculation
    print("This script calculates the expected chain quality Q based on the formula:")
    print("Q = (α + 2α²β - α²βp) / (αβ + α + β/α)\n")
    print(f"Given values:")
    print(f"  Adversary mining power (β) = {beta}")
    print(f"  Honest mining power (α) = 1 - {beta} = {alpha}")
    print(f"  Tie-breaking probability (p) = {p}\n")
    print("Substituting these values into the formula:\n")

    # Print the full equation with substituted values
    num_str = f"({alpha:.3f} + 2*({alpha:.3f}**2)*{beta:.3f} - ({alpha:.3f}**2)*{beta:.3f}*{p:.3f})"
    den_str = f"({alpha:.3f}*{beta:.3f} + {alpha:.3f} + {beta:.3f}/{alpha:.3f})"
    
    num_val_part1 = alpha
    num_val_part2 = 2 * (alpha**2) * beta
    num_val_part3 = (alpha**2) * beta * p
    den_val_part1 = alpha * beta
    den_val_part2 = alpha
    den_val_part3 = beta / alpha
    
    print(f"Q = {num_str} / {den_str}")
    print(f"Q = ({num_val_part1:.4f} + {num_val_part2:.4f} - {num_val_part3:.4f}) / ({den_val_part1:.4f} + {den_val_part2:.4f} + {den_val_part3:.4f})")
    print(f"Q = {numerator:.4f} / {denominator:.4f}\n")


    # Print the final result
    print(f"The expected chain quality is: {chain_quality:.4f}")


# --- Main execution ---
# You can change these values to explore different scenarios
adversary_power_beta = 0.3  # Example: Adversary has 30% of the mining power
tie_break_prob_p = 0.5      # Example: 50% chance to favor the adversary in a tie

calculate_chain_quality(adversary_power_beta, tie_break_prob_p)