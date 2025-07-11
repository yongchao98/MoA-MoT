import sys

def solve():
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    The script takes two command-line arguments:
    1. beta (β): The fraction of mining power controlled by the adversary (0 < beta < 0.5).
    2. p: The probability that honest miners choose the adversary's block in a tie (0 <= p <= 1).
    """
    # Check for correct number of command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <beta> <p>")
        print("  <beta>: adversary's mining power, a value between 0 and 0.5.")
        print("  <p>: tie-breaking probability, a value between 0 and 1.")
        return

    # Parse and validate inputs
    try:
        beta = float(sys.argv[1])
        p = float(sys.argv[2])
    except ValueError:
        print("Error: Invalid input. beta and p must be numeric values.")
        return

    # The underlying Markov model requires β < 0.5 for the state probabilities to converge.
    if not (0 < beta < 0.5):
        print(f"Error: The input for beta (β) must be strictly between 0 and 0.5 for this model.")
        if beta <= 0:
            quality = 1.0
            print("With beta <= 0, there is no adversary, so the chain quality is 1.")
            print(f"<<<{quality}>>>")
        else: # beta >= 0.5
            quality = 0.0
            print("With beta >= 0.5, the adversary's private chain grows indefinitely long.")
            print("The longest chain contains no honest blocks, so the quality is 0.")
            print(f"<<<{quality}>>>")
        return
        
    if not (0 <= p <= 1):
        print("Error: The input for p must be between 0 and 1 (inclusive).")
        return

    # Let alpha be the honest miners' power portion
    alpha = 1.0 - beta

    # The formula for expected chain quality is derived from a Markov chain model of the system.
    # It represents the long-term ratio of the rate of honest blocks being added to the chain
    # to the rate of total blocks being added.
    # Formula: Quality = ( (1 - β) * (1 + 3β - pβ) ) / (1 + β)^2
    
    numerator_val = alpha * (1 + 3*beta - p*beta)
    denominator_val = (1 + beta)**2
    
    chain_quality = numerator_val / denominator_val

    # Output the explanation and result
    print("This script calculates the expected chain quality based on the provided model.")
    print(f"Inputs: Adversary mining power (beta) = {beta}, Tie-breaking probability (p) = {p}")
    print("-" * 20)

    print("The derived formula for chain quality is:")
    print("Quality = ( (1 - beta) * (1 + 3*beta - p*beta) ) / ( (1 + beta)^2 )")
    print("\nPlugging in the given values:")

    print("\n1. Calculate the numerator:")
    print(f"  (1 - {beta}) * (1 + 3*{beta} - {p}*{beta})")
    print(f"= {alpha} * ({1 + 3*beta - p*beta})")
    print(f"= {numerator_val}")
    
    print("\n2. Calculate the denominator:")
    print(f"  (1 + {beta})^2")
    print(f"= ({1+beta})^2")
    print(f"= {denominator_val}")
    
    print("\n3. Calculate the final Chain Quality:")
    print(f"  Quality = Numerator / Denominator")
    print(f"  Quality = {numerator_val} / {denominator_val}")
    print(f"  Quality = {chain_quality}")
    
    print("\n" + "="*40)
    print("Final Result")
    print(f"The expected chain quality is: {chain_quality}")
    print("="*40)
    
    # Final answer in the specified format
    print(f"<<<{chain_quality}>>>")


if __name__ == '__main__':
    solve()
