import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin-like system
    under a selfish mining attack with specific tie-breaking rules.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability that honest miners choose the adversary's block
                     in a tie-breaking situation (0 <= p <= 1).

    Returns:
        float: The expected chain quality, or None if inputs are invalid.
    """
    # Validate inputs
    if not (0 < beta < 1):
        print("Error: Adversary's mining power (beta) must be between 0 and 1.", file=sys.stderr)
        return None
    if not (0 <= p <= 1):
        print("Error: Tie-breaking probability (p) must be between 0 and 1.", file=sys.stderr)
        return None

    # The formula for chain quality is derived from a Markov chain model of the selfish mining process.
    # It is the ratio of the expected number of honest blocks to the expected total number of blocks
    # added to the longest chain over a long period.
    
    # Numerator of the chain quality formula. This is proportional to the
    # long-term expected number of honest blocks added to the chain.
    # Formula: (1-β)² * (1 + 2β - 2β² - p*β*(1-β))
    honest_blocks_numerator = (1 - beta)**2 * (1 + 2 * beta - 2 * beta**2 - p * beta * (1 - beta))

    # Denominator of the chain quality formula. This is proportional to the
    # long-term expected length of the chain (total blocks).
    # Formula: 1 - β² + β³
    total_blocks_denominator = 1 - beta**2 + beta**3

    # The chain quality is the ratio of the two values.
    chain_quality = honest_blocks_numerator / total_blocks_denominator

    print("Inputs:")
    print(f"  Adversary's mining power (β): {beta}")
    print(f"  Tie-breaking probability (p): {p}\n")
    
    print("The final equation for Chain Quality is Numerator / Denominator.")
    print("---------------------------------------------------------------")
    print(f"Numerator (proportional to honest blocks): (1-{beta})² * (1 + 2*{beta} - 2*({beta}²) - {p}*{beta}*(1-{beta})) = {honest_blocks_numerator}")
    print(f"Denominator (proportional to total blocks): 1 - ({beta}²) + ({beta}³) = {total_blocks_denominator}")
    print("---------------------------------------------------------------")
    print(f"Expected Chain Quality = {honest_blocks_numerator} / {total_blocks_denominator}")
    print(f"                     = {chain_quality}")
    
    return chain_quality


if __name__ == '__main__':
    # --- User-configurable values ---
    # β (beta): The fraction of the total mining power controlled by the adversary.
    # For selfish mining to be effective, beta typically needs to be > 1/4. Let's use 1/3 as an example.
    adversary_power_beta = 1/3

    # p: The probability that honest miners choose the adversary's block to extend
    # when a tie occurs. A value of 0.5 means they choose randomly.
    tie_break_prob_p = 0.5
    
    calculate_chain_quality(adversary_power_beta, tie_break_prob_p)
    # The final answer will be printed by the function.
    # We can capture it if needed.
    result = calculate_chain_quality(adversary_power_beta, tie_break_prob_p)
    if result is not None:
        # The instruction asks for the answer in a specific format at the end.
        final_answer = round(result, 6) # Round for a clean final output
        # Let's hide the final answer wrapper as the instructions state not to include multiple code blocks.
        # But for the purpose of the exercise, the numerical value is computed and can be extracted.
