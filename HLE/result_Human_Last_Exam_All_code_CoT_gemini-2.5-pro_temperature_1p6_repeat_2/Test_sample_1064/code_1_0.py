import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a selfish mining attack.

    Args:
        beta (float): The proportion of mining power held by the adversary (0 <= beta < 1).
        p (float): The probability an honest miner chooses the adversary's block in a tie (0 <= p <= 1).
    """
    # Check for valid inputs
    if not (0 <= beta < 1):
        print("Error: beta (adversary's mining power) must be between 0 (inclusive) and 1 (exclusive).", file=sys.stderr)
        return
    if not (0 <= p <= 1):
        print("Error: p (tie-breaking probability) must be between 0 and 1 (inclusive).", file=sys.stderr)
        return

    # The derived formula for Chain Quality (CQ) is:
    # CQ = (Numerator) / (Denominator)
    # Numerator = (1 - beta)^2 * (1 + (2-p)*beta - (2-p)*beta^2)
    # Denominator = 1 - beta^2 + beta^3

    # --- Calculation ---

    # Calculate numerator components
    term1_num_val = (1 - beta)**2
    two_minus_p = 2 - p
    term2_num_val = 1 + two_minus_p * beta - two_minus_p * (beta**2)
    numerator = term1_num_val * term2_num_val

    # Calculate denominator
    denominator = 1 - beta**2 + beta**3

    # Calculate final chain quality
    # Denominator is guaranteed to be non-zero for 0 <= beta < 1
    chain_quality = numerator / denominator

    # --- Output ---

    print(f"For an adversary with mining power beta = {beta} and a tie-breaking probability p = {p}:\n")
    print("The expected chain quality is calculated using the formula:")
    print("CQ = ((1 - beta)^2 * (1 + (2-p)*beta - (2-p)*beta^2)) / (1 - beta^2 + beta^3)\n")
    
    print("Step-by-step calculation:")
    
    # Numerator calculation breakdown
    print(f"Numerator = (1 - {beta})^2 * (1 + (2-{p})*{beta} - (2-{p})*({beta}^2))")
    print(f"Numerator = ({1-beta})^2 * (1 + {two_minus_p}*{beta} - {two_minus_p}*{beta**2})")
    print(f"Numerator = {term1_num_val} * ({term2_num_val})")
    print(f"Numerator = {numerator}\n")
    
    # Denominator calculation breakdown
    print(f"Denominator = 1 - {beta}^2 + {beta}^3")
    print(f"Denominator = 1 - {beta**2} + {beta**3}")
    print(f"Denominator = {denominator}\n")
    
    # Final Result
    print(f"Expected Chain Quality = Numerator / Denominator")
    print(f"Expected Chain Quality = {numerator} / {denominator}")
    
    # Use f-string formatting to control the number of decimal places for the final answer
    final_answer = f"{chain_quality:.6f}"
    print(f"Expected Chain Quality = {final_answer}")
    
    # Return final answer in the requested format
    return final_answer


if __name__ == '__main__':
    # --- Parameters ---
    # Adversary's mining power proportion.
    # Selfish mining is generally profitable for beta > 1/3 (depending on p).
    BETA = 0.35
    
    # Probability that honest miners choose the adversary's block when there is a tie.
    # p=0.5 represents random choice.
    P_TIE = 0.5
    
    result = calculate_chain_quality(BETA, P_TIE)
    # The final answer is wrapped for easy parsing.
    if result is not None:
        print(f"\n<<<{result}>>>")