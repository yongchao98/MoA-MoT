import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a specific selfish mining scenario.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability an honest miner chooses the adversary's block in a tie (0 <= p <= 1).
    """
    if not (0 < beta < 1):
        print("Error: Adversary's mining power (beta) must be between 0 and 1.", file=sys.stderr)
        return
    if not (0 <= p <= 1):
        print("Error: The probability (p) must be between 0 and 1.", file=sys.stderr)
        return

    # Calculate the expected number of honest blocks added to the chain per cycle.
    # Formula: E[h] = (1 - beta) * (1 + beta * (1 - beta) * (2 - p))
    e_h_term1 = 1 - beta
    e_h_term2 = 1 + beta * (1 - beta) * (2 - p)
    e_h = e_h_term1 * e_h_term2

    # Calculate the expected total number of blocks added to the chain per cycle.
    # Formula: E[L] = (1 - beta^2 + beta^3) / (1 - beta)
    e_l_numerator = 1 - beta**2 + beta**3
    e_l_denominator = 1 - beta
    e_l = e_l_numerator / e_l_denominator

    # Chain quality is the ratio of expected honest blocks to expected total blocks.
    chain_quality = e_h / e_l
    
    # The final equation is Q = E[h] / E[L]
    # We can also write it as Q = ((1-beta)^2 * (1 + beta*(1-beta)*(2-p))) / (1 - beta^2 + beta^3)
    final_numerator = (1 - beta)**2 * (1 + beta * (1 - beta) * (2 - p))
    final_denominator = 1 - beta**2 + beta**3

    print(f"For beta = {beta} and p = {p}:")
    print("\n--- Calculation Steps ---")
    print(f"Numerator of the chain quality formula (E[h] * (1-beta)): {final_numerator}")
    print(f"Denominator of the chain quality formula (E[L] * (1-beta)): {final_denominator}")
    print("\nFinal Equation:")
    print(f"Chain Quality = {final_numerator} / {final_denominator}")
    print(f"\nExpected Chain Quality: {chain_quality}")
    
    # Printing the final answer in the required format
    print(f"\n<<<{chain_quality}>>>")

if __name__ == '__main__':
    # --- Parameters ---
    # The portion of mining power controlled by the selfish miner.
    # Must be between 0 and 1 (exclusive).
    adversary_power_beta = 0.3

    # The probability that an honest miner will choose to mine on the
    # adversary's block in case of a tie.
    # Must be between 0 and 1 (inclusive).
    tie_preference_p = 0.5

    calculate_chain_quality(adversary_power_beta, tie_preference_p)