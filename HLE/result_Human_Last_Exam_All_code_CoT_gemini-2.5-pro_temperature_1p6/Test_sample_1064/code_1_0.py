import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin-like system
    with a selfish miner.

    Args:
        beta (float): The portion of mining power held by the selfish miner (0 < beta < 1).
        p (float): The probability an honest miner chooses the attacker's block in a tie (0 <= p <= 1).
        
    Returns:
        float: The expected chain quality.
    """
    if not (0 < beta < 1):
        raise ValueError("beta must be between 0 and 1.")
    if not (0 <= p <= 1):
        raise ValueError("p must be between 0 and 1.")

    # Let alpha = 1 - beta
    alpha = 1 - beta

    # Numerator of the chain quality formula
    # (1-beta)^2 * [1 + beta*(1-beta)*(2-p)]
    numerator = (alpha**2) * (1 + beta * alpha * (2 - p))

    # Denominator of the chain quality formula
    # 1 - beta^2 + beta^3
    denominator = 1 - beta**2 + beta**3

    if denominator == 0:
        return float('inf') # Should not happen for 0 < beta < 1
        
    chain_quality = numerator / denominator
    return chain_quality, numerator, denominator

# --- Parameters for the calculation ---
# beta: adversary's mining power portion.
# p: probability honest miners choose adversary's block in a tie.
# We assume beta=0.4 and p=0.5 as illustrative values.
beta = 0.4
p = 0.5

# Calculate the chain quality and the components of the formula
try:
    quality, num_val, den_val = calculate_chain_quality(beta, p)

    # --- Output the results as requested ---
    print("Calculating the expected chain quality with the following parameters:")
    print(f"Adversary's mining power (β): {beta}")
    print(f"Tie-breaking probability (p): {p}")
    print("-" * 30)

    # Print the final equation with all the numbers
    print("Final Equation:")
    # The double curly braces {{ and }} are used to print literal curly braces in an f-string
    print(f"Chain Quality = ((1 - {beta})^2 * (1 + {beta} * (1 - {beta}) * (2 - {p}))) / (1 - {beta}^2 + {beta}^3)")
    
    # Print the evaluated numerator and denominator
    print("\nCalculation Steps:")
    print(f"Numerator   = ({1-beta:.2f}^2) * (1 + {beta} * {1-beta:.2f} * {2-p:.2f}) = {num_val}")
    print(f"Denominator = 1 - {beta**2:.3f} + {beta**3:.3f} = {den_val}")
    
    # Print the final result
    print("\nResult:")
    print(f"Expected Chain Quality = {num_val} / {den_val} = {quality}")

except ValueError as e:
    print(f"Error: {e}", file=sys.stderr)

# The final answer in the required format
# print(f"\n<<<{quality}>>>")