import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability an honest miner chooses the adversary's block in a tie (0 <= p <= 1).
    """

    if not (0 < beta < 1):
        print("Error: beta (β) must be between 0 and 1.", file=sys.stderr)
        return
    if not (0 <= p <= 1):
        print("Error: p must be between 0 and 1.", file=sys.stderr)
        return

    print(f"Calculating chain quality for β = {beta} and p = {p}")
    print("-" * 30)

    # If the adversary's power is 50% or more, their chain will grow indefinitely.
    # In the long run, no honest blocks will be in the main chain.
    if beta >= 0.5:
        print("Since β >= 0.5, the adversary will dominate the chain.")
        print("Expected Chain Quality = 0")
        print("\nFinal Answer:")
        print("<<<0>>>")
        return

    # Numerator components
    num1 = 1 - 2 * beta
    num2 = 1 - beta
    num3 = 1 + 2 * beta - p * beta
    
    # Denominator
    den = 1 - beta - beta**2

    if den == 0:
        print("Error: Denominator is zero. Cannot calculate for this beta value.", file=sys.stderr)
        return
        
    numerator_val = num1 * num2 * num3
    chain_quality = numerator_val / den

    print("The formula for chain quality (CQ) is:")
    print("CQ = ((1 - 2β) * (1 - β) * (1 + 2β - pβ)) / (1 - β - β²)\n")
    
    print("Plugging in the values:")
    print(f"Numerator = ({num1:.4f}) * ({num2:.4f}) * ({num3:.4f}) = {numerator_val:.4f}")
    print(f"Denominator = {den:.4f}")
    
    print(f"\nFinal Equation:")
    print(f"Chain Quality = ({num1:.4f} * {num2:.4f} * {num3:.4f}) / {den:.4f} = {chain_quality}")

    # For final answer submission format
    final_answer_str = f"<<<{chain_quality}>>>"
    print(f"\nFinal Answer in submission format:\n{final_answer_str}")


if __name__ == '__main__':
    # Example values for demonstration
    # You can change these values to test other scenarios
    beta_value = 0.3
    p_value = 0.5
    
    calculate_chain_quality(beta_value, p_value)
