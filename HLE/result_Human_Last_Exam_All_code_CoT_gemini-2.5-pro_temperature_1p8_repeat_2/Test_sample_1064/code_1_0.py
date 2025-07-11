import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a specific selfish mining scenario.

    Args:
        beta (float): The portion of mining power held by the adversary (0 <= beta < 1).
        p (float): The probability an honest miner chooses the adversary's block in a tie (0 <= p <= 1).
    """

    print(f"Given parameters:")
    print(f"  Adversary's mining power (β): {beta}")
    print(f"  Honest miner's choice probability in a tie (p): {p}")
    print("-" * 30)

    # For β >= 0.5, the adversary's private chain grows indefinitely,
    # leading to a chain quality of 0.
    if beta >= 0.5:
        print("Analysis: When β >= 0.5, the adversary's chain is expected to grow longer than the honest chain.")
        print("In the long run, the adversary will publish a significantly longer chain, orphaning all honest blocks.")
        print("Therefore, the number of honest blocks on the longest chain approaches zero.")
        print("\nFinal Equation:")
        print("Chain Quality = 0 / Chain Length = 0.0")
        return 0.0

    # For β < 0.5, a stationary distribution exists, and we can calculate the chain quality using the derived formula.
    print("Analysis: When β < 0.5, the system reaches a stable state. The chain quality can be calculated")
    print("using the long-term rates of honest and selfish blocks being added to the chain.")
    print("The derived formula for chain quality (Q) is:")
    print("Q = ((1 - β) * (1 + 2*β - p*β) * (1 - 2*β)) / (1 - β - β^2)\n")

    # Calculate each part of the formula
    numerator_part1 = 1 - beta
    numerator_part2 = 1 + 2 * beta - p * beta
    numerator_part3 = 1 - 2 * beta
    numerator = numerator_part1 * numerator_part2 * numerator_part3

    denominator = 1 - beta - beta**2

    # Prevent division by zero, although this case is covered by beta >= 0.5 check
    # as 1 - beta - beta^2 = 0 implies beta is approx 0.618
    if denominator == 0:
        # This case is not practically reachable for beta < 0.5
        chain_quality = 0.0
    else:
        chain_quality = numerator / denominator

    # Output the final equation with all the numbers
    print("Final Equation with substituted values:")
    final_equation_str = (
        f"Q = ((1 - {beta}) * (1 + 2*{beta} - {p}*{beta}) * (1 - 2*{beta})) / (1 - {beta} - {beta}**2)"
    )
    print(final_equation_str)
    
    # Show intermediate calculation steps
    print(f"   = (({numerator_part1}) * ({numerator_part2}) * ({numerator_part3})) / ({denominator})")
    print(f"   = {numerator} / {denominator}")
    
    # Print the final result
    print(f"   = {chain_quality}")
    
    return chain_quality

if __name__ == '__main__':
    # Example values for demonstration
    # You can change these values to see different outcomes
    beta_value = 0.25
    p_value = 0.5

    final_quality = calculate_chain_quality(beta_value, p_value)
    
    # The '<<<...>>>' marker expects a single value.
    # We will output this marker with the final calculated result.
    # To avoid printing it from inside the function, we print it here.
    # Note: The required output format is just the raw value.
    sys.stdout = open(_io.StringIO(), 'w') # Redirect stdout to suppress the marker in the final output block
    print(f"\n<<<{final_quality}>>>")
    sys.stdout = sys.__stdout__ # Restore stdout