import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a specific selfish mining model.

    Args:
        beta (float): The adversary's portion of the total mining power. Must be between 0 and 1.
        p (float): The probability that an honest miner chooses the adversary's block in a tie. Must be between 0 and 1.
    """
    if not (0 <= beta < 1):
        print("Error: beta (β) must be between 0 and 1 (exclusive of 1).")
        return
    if not (0 <= p <= 1):
        print("Error: p must be between 0 and 1.")
        return

    print(f"Calculating expected chain quality for β = {beta} and p = {p}\n")

    # If the adversary has half or more of the mining power, their lead will
    # grow indefinitely in the long run. They can orphan all honest blocks.
    if beta >= 0.5:
        print("When β >= 0.5, the adversary's lead grows indefinitely.")
        print("The expected chain quality approaches 0.")
        print("\nFinal Answer: Expected Chain Quality = 0")
        return

    # For β < 0.5, we use the derived formula:
    # q = ( (1-β)(1-2β)(1 + 2β - βp) ) / ( 1 - β - β^2 )

    # Calculate each term of the formula to show the steps
    term1 = 1 - beta
    term2 = 1 - 2 * beta
    term3 = 1 + 2 * beta - beta * p
    
    numerator = term1 * term2 * term3
    denominator = 1 - beta - beta**2
    
    if denominator == 0:
        # This case occurs at β = (sqrt(5)-1)/2 ≈ 0.618, which is already
        # handled by the β >= 0.5 check. This is for numerical safety.
        print("Error: Denominator is zero. This case should not be reached for β < 0.5.")
        return
        
    quality = numerator / denominator

    print("The formula for expected chain quality (q) is:")
    print("q = [ (1-β) * (1-2β) * (1 + 2β - βp) ] / [ 1 - β - β^2 ]\n")

    print("Substituting the given values:")
    print(f"q = [ ({term1:.4f}) * ({term2:.4f}) * ({term3:.4f}) ] / [ {denominator:.4f} ]")
    print(f"q = [ {numerator:.4f} ] / [ {denominator:.4f} ]\n")

    print(f"Final Answer: Expected Chain Quality = {quality:.6f}")


if __name__ == '__main__':
    # Example usage. You can change these values.
    # Default values if not provided via command line
    beta_val = 0.333
    p_val = 0.5

    if len(sys.argv) == 3:
        try:
            beta_val = float(sys.argv[1])
            p_val = float(sys.argv[2])
        except ValueError:
            print("Invalid input. Please provide two numbers for beta and p.")
            print("Usage: python your_script_name.py <beta> <p>")
            sys.exit(1)
            
    calculate_chain_quality(beta_val, p_val)
