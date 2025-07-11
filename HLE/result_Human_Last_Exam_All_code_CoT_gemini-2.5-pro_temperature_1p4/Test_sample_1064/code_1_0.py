import sys

def calculate_chain_quality():
    """
    Calculates the expected chain quality in Bitcoin under a selfish mining attack.
    Prompts the user for the adversary's mining power (beta) and the tie-breaking probability (p).
    """
    try:
        beta_str = input("Enter the adversary's mining power proportion (β, e.g., 0.3): ")
        beta = float(beta_str)
        p_str = input("Enter the probability honest miners choose the adversary's block in a tie (p, e.g., 0.5): ")
        p = float(p_str)
    except ValueError:
        print("Invalid input. Please enter numerical values.", file=sys.stderr)
        sys.exit(1)

    if not (0 <= beta < 1):
        print("Error: Beta (β) must be a value between 0 and 1.", file=sys.stderr)
        sys.exit(1)
    if not (0 <= p <= 1):
        print("Error: Probability (p) must be a value between 0 and 1.", file=sys.stderr)
        sys.exit(1)

    # If beta is 0.5 or greater, the selfish miner's chain grows indefinitely.
    # The long-run chain quality is 0.
    if beta >= 0.5:
        chain_quality = 0.0
        print("\nFor β >= 0.5, the adversary will always have a longer chain in the long run.")
        print("Expected Chain Quality = 0.0")
        print(f"\nFinal Answer: {chain_quality}")
        print(f'<<<{chain_quality}>>>')
        return

    # For beta < 0.5, we use the derived formula.
    # Numerator calculation
    term1_num = 1 - beta
    term2_num = 1 - (2 * beta)
    term3_num = 1 + (2 * beta) - (beta * p)
    numerator = term1_num * term2_num * term3_num

    # Denominator calculation
    denominator = 1 - beta - (beta**2)

    # Denominator can't be zero for beta < 0.5
    chain_quality = numerator / denominator

    print("\nBased on the analysis of the selfish mining stochastic process, the expected chain quality is calculated as:")
    print("CQ = ((1 - β) * (1 - 2*β) * (1 + 2*β - β*p)) / (1 - β - β²)")
    
    print("\nPlugging in the values:")
    # Fulfills the requirement to "output each number in the final equation"
    print(f"Numerator = (1 - {beta}) * (1 - 2*{beta}) * (1 + 2*{beta} - {beta}*{p})")
    print(f"Numerator = ({term1_num}) * ({term2_num}) * ({term3_num}) = {numerator}")
    print(f"\nDenominator = 1 - {beta} - {beta}²")
    print(f"Denominator = {denominator}")
    
    print(f"\nCQ = {numerator} / {denominator}")
    
    print(f"\nFinal Answer: {chain_quality}")
    print(f'<<<{chain_quality}>>>')

if __name__ == '__main__':
    calculate_chain_quality()