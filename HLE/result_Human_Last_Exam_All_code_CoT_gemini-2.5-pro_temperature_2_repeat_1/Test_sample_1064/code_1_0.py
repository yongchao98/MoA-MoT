import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality based on the provided model.

    Args:
        beta (float): The adversary's proportion of the total mining power (0 < beta < 0.5).
        p (float): The probability that honest miners choose the adversary's block in a tie (0 <= p <= 1).
    """
    # Parameter validation
    if not (0 < beta < 1):
        print("Error: beta must be a value between 0 and 1 (exclusive).", file=sys.stderr)
        return
        
    if not (0 <= p <= 1):
        print("Error: p must be a value between 0 and 1 (inclusive).", file=sys.stderr)
        return
    
    if beta >= 0.5:
        print("\nWarning: For beta >= 0.5, the selfish mining strategy is always dominant.", file=sys.stderr)
        print("The chain quality will approach 0 as the attacker can always build a longer chain.", file=sys.stderr)

    # --- Step 1: Calculate the expected number of honest blocks per cycle (E_h) ---
    # Formula: E_h = (1-β) * (1 + β(1-β)(2-p))
    term_1_beta = 1 - beta
    term_beta_p = beta * term_1_beta * (2 - p)
    expected_honest_blocks = term_1_beta * (1 + term_beta_p)

    # --- Step 2: Calculate the expected total blocks added to the chain per cycle (E_L) ---
    # This is the sum of expected honest blocks and expected adversary blocks.
    # Formula: E_L = (1 - β^2 + β^3) / (1-β)
    numerator = 1 - beta**2 + beta**3
    denominator = 1 - beta
    if denominator == 0:
        total_expected_blocks = float('inf')
    else:
        total_expected_blocks = numerator / denominator

    # --- Step 3: Calculate the final Chain Quality (Q) ---
    # Formula: Q = E_h / E_L
    if total_expected_blocks == 0 or total_expected_blocks == float('inf'):
        chain_quality = 0.0
    else:
        chain_quality = expected_honest_blocks / total_expected_blocks

    # --- Step 4: Output the results as an equation ---
    print("\n--- Calculation Results ---")
    print(f"Given beta = {beta} and p = {p}:")
    
    print("\nThe final equation for Chain Quality is:")
    print("Chain Quality = (Expected Number of Honest Blocks) / (Expected Total Chain Length)")
    print(f"Chain Quality = {expected_honest_blocks} / {total_expected_blocks}")
    print(f"Chain Quality = {chain_quality}")


if __name__ == '__main__':
    # You can change these values to test different scenarios
    # Adversary's mining power (beta)
    adversary_power_beta = 0.3
    # Tie-breaking probability (p)
    tie_breaking_p = 0.5
    
    calculate_chain_quality(adversary_power_beta, tie_breaking_p)
    
    # The final answer format is just the numerical value
    # Let's re-run to capture it for the final output tag
    term_1_beta = 1 - adversary_power_beta
    term_beta_p = adversary_power_beta * term_1_beta * (2 - tie_breaking_p)
    expected_honest_blocks = term_1_beta * (1 + term_beta_p)
    numerator = 1 - adversary_power_beta**2 + adversary_power_beta**3
    denominator = 1 - adversary_power_beta
    total_expected_blocks = numerator / denominator
    chain_quality = expected_honest_blocks / total_expected_blocks
    print(f"\n<<<{chain_quality}>>>")
