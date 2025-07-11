import sys

def solve_and_print_chain_quality(beta, p):
    """
    Calculates and prints the expected chain quality for the described selfish mining scenario.
    
    Args:
        beta (float): The portion of mining power of the adversary (0 < beta < 1).
        p (float): The probability of honest miners choosing the adversary's block in a tie (0 <= p <= 1).
        
    Returns:
        float: The calculated expected chain quality.
    """
    print("--- Chain Quality Calculation ---")
    print(f"Inputs: Adversary Power (β) = {beta:.4f}, Tie-breaking prob (p) = {p:.4f}\n")

    # If the adversary's mining power is 0.5 or greater, their private chain is expected
    # to grow indefinitely, making the chain quality (honest blocks / total blocks) approach 0.
    if beta >= 0.5:
        print("Result: For β >= 0.5, the selfish miner's chain grows indefinitely longer.")
        print("Expected Chain Quality = 0.0")
        return 0.0

    # For β < 0.5, the chain quality (CQ) is calculated using the derived formula.
    print("The formula for Chain Quality (CQ) for β < 0.5 is:")
    print("CQ = [ (1 - β) * (1 + 3*β - p*β) ] / (1 + β)²\n")
    
    alpha = 1.0 - beta
    
    # Numerator calculation
    numerator_term1 = alpha
    numerator_term2 = 1 + 3 * beta - p * beta
    numerator = numerator_term1 * numerator_term2
    
    # Denominator calculation
    denominator_base = 1 + beta
    denominator = denominator_base ** 2
    
    chain_quality = numerator / denominator
    
    # As requested, we print each number in the final equation.
    print("Calculation Breakdown:")
    print(f"  Numerator = (1 - {beta:.4f}) * (1 + 3*{beta:.4f} - {p:.4f}*{beta:.4f})")
    print(f"            = ({numerator_term1:.4f}) * ({numerator_term2:.4f})")
    print(f"            = {numerator:.4f}")
    
    print(f"\n  Denominator = (1 + {beta:.4f})²")
    print(f"              = ({denominator_base:.4f})²")
    print(f"              = {denominator:.4f}")
    
    print(f"\nFinal Result: CQ = {numerator:.4f} / {denominator:.4f}")
    print(f"Expected Chain Quality = {chain_quality:.4f}")
    
    return chain_quality

# --- User-configurable values ---
# You can change these values to see the result for different scenarios.

# beta: The portion of mining power of the adversary (must be between 0 and 1).
# Using the example value 1/3.
beta = 1/3

# p: The probability of honest miners choosing the adversary's block in a tie.
# Using the example value 0.5.
p = 0.5
# ---------------------------------

# Basic input validation
if not (0 < beta < 1):
    print("Error: beta must be between 0 and 1.", file=sys.stderr)
elif not (0 <= p <= 1):
    print("Error: p must be between 0 and 1.", file=sys.stderr)
else:
    solve_and_print_chain_quality(beta, p)