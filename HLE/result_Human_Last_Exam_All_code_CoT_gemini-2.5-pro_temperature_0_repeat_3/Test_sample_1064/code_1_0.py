import math

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power held by the adversary (0 <= beta < 1).
        p (float): The probability that honest miners choose the adversary's block in a tie (0 <= p <= 1).
    """
    print(f"Calculating chain quality for β = {beta} and p = {p}\n")

    # If the adversary's power is 0.5 or greater, their lead grows indefinitely.
    # Honest blocks are never part of the longest chain in the long run.
    if beta >= 0.5:
        print("Since β is 0.5 or greater, the adversary's chain will always be longer.")
        print("The expected chain quality is 0.")
        final_quality = 0
        print("\n--- Final Result ---")
        print(f"Expected Chain Quality: {final_quality}")
        return final_quality

    # The formula is derived from a Markov chain model of selfish mining.
    # It is valid for 0 <= beta < 0.5.
    print("The formula for chain quality (q) when β < 0.5 is:")
    print("q = ((1 - β) * (1 + 2*β - p*β) * (1 - 2*β)) / (1 - β - β^2)\n")

    # Calculate numerator
    term1_num = 1 - beta
    term2_num = 1 + 2 * beta - p * beta
    term3_num = 1 - 2 * beta
    numerator = term1_num * term2_num * term3_num

    # Calculate denominator
    denominator = 1 - beta - beta**2
    
    if denominator <= 0:
        print("Error: The denominator is zero or negative, which should not happen for β < 0.5.")
        print("Please check the value of β. It must be less than (sqrt(5)-1)/2 ≈ 0.618 for a meaningful result in this model, and less than 0.5 for stability.")
        return None

    # Calculate final quality
    quality = numerator / denominator

    # Print the step-by-step calculation
    print("--- Calculation Steps ---")
    print("1. Calculate the numerator:")
    print(f"   Numerator = (1 - {beta}) * (1 + 2*{beta} - {p}*{beta}) * (1 - 2*{beta})")
    print(f"   Numerator = ({term1_num}) * ({term2_num}) * ({term3_num})")
    print(f"   Numerator = {numerator:.4f}")

    print("\n2. Calculate the denominator:")
    print(f"   Denominator = 1 - {beta} - {beta}^2")
    print(f"   Denominator = 1 - {beta} - {beta**2:.4f}")
    print(f"   Denominator = {denominator:.4f}")

    print("\n3. Calculate the final chain quality:")
    print(f"   Chain Quality = Numerator / Denominator")
    print(f"   Chain Quality = {numerator:.4f} / {denominator:.4f}")
    
    print("\n--- Final Result ---")
    print(f"Expected Chain Quality: {quality:.4f}")
    
    return quality

if __name__ == '__main__':
    # --- User-defined parameters ---
    # β (beta): Adversary's portion of mining power. Must be between 0 and 1.
    # The model is stable for beta < 0.5.
    adversary_power_beta = 0.3

    # p: Probability honest miners choose the adversary's block in a tie.
    # Must be between 0 and 1.
    tie_breaking_prob_p = 0.5
    # -----------------------------

    final_answer = calculate_chain_quality(adversary_power_beta, tie_breaking_prob_p)
    
    # The following line is for the final answer submission format.
    if final_answer is not None:
        print(f"\n<<<__{final_answer:.4f}__>>>")
