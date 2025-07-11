import sys

def calculate_chain_quality():
    """
    Calculates the expected chain quality under a specific selfish mining model.

    The user is prompted to enter the adversary's mining power (beta) and
    the probability of choosing the adversary's block in a tie (p).

    The function then calculates and prints the expected chain quality based on
    a pre-derived formula from a Markov chain model of the system.
    """
    try:
        beta_str = input("Enter the adversary's mining power fraction (beta, e.g., 0.3): ")
        beta = float(beta_str)
        if not (0 <= beta < 1):
            print("Error: beta must be between 0 and 1.", file=sys.stderr)
            return

        p_str = input("Enter the probability of choosing adversary's block in a tie (p, e.g., 0.5): ")
        p = float(p_str)
        if not (0 <= p <= 1):
            print("Error: p must be between 0 and 1.", file=sys.stderr)
            return

    except ValueError:
        print("Error: Invalid input. Please enter numerical values.", file=sys.stderr)
        return

    print("\n--- Calculation ---")
    
    # If beta >= 0.5, the selfish miner's lead grows indefinitely.
    # The Markov chain model does not converge, and they take over the chain.
    if beta >= 0.5:
        print(f"For beta = {beta}, the adversary's power is too high.")
        print("The expected chain quality approaches 0.")
        chain_quality = 0
    else:
        # These are the components of the derived formula:
        # Numerator = (1-beta) * (1 + 2*beta - p*beta) * (1-2*beta)
        # Denominator = 1 - beta - beta**2

        term1 = 1 - beta
        term2 = 1 + 2 * beta - p * beta
        term3 = 1 - 2 * beta
        numerator = term1 * term2 * term3

        denominator = 1 - beta - beta**2
        
        # Output the equation with the user-provided numbers
        print("The formula for Chain Quality (CQ) is:")
        print("CQ = ((1 - beta) * (1 + 2*beta - p*beta) * (1 - 2*beta)) / (1 - beta - beta^2)\n")

        print("Substituting the values:")
        print(f"CQ = ((1 - {beta}) * (1 + 2*{beta} - {p}*{beta}) * (1 - 2*{beta})) / (1 - {beta} - {beta}**2)")
        
        print("\nStep-by-step evaluation:")
        print(f"Term 1 (1 - beta): {term1}")
        print(f"Term 2 (1 + 2*beta - p*beta): {term2}")
        print(f"Term 3 (1 - 2*beta): {term3}")
        print(f"Numerator (Term 1 * Term 2 * Term 3): {numerator}")
        print(f"Denominator (1 - beta - beta^2): {denominator}\n")
        
        if denominator == 0:
             print("Error: Calculation resulted in division by zero.", file=sys.stderr)
             return
             
        chain_quality = numerator / denominator

    print("---------------------------------")
    print(f"Expected Chain Quality: {chain_quality}")
    print("---------------------------------")


if __name__ == "__main__":
    calculate_chain_quality()