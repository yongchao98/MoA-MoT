def calculate_chain_quality():
    """
    Calculates the expected chain quality in a Bitcoin system
    under a specific selfish mining attack scenario.

    The user can modify the hardcoded values for beta and p below.
    """
    # Adversary's portion of mining power.
    beta = 0.3
    # Probability that honest miners choose the adversary's block in a tie.
    p = 0.5

    print(f"Analyzing the system with the following parameters:")
    print(f"Adversary's mining power (β): {beta}")
    print(f"Honest miners' tie-breaking preference (p): {p}\n")

    # Validate input parameters
    if not (0 <= beta < 1):
        print("Error: Adversary mining power (β) must be a value between 0 and 1.")
        return
    if not (0 <= p <= 1):
        print("Error: Tie-breaking probability (p) must be a value between 0 and 1.")
        return

    # If the adversary's power is 50% or more, they can always ensure their
    # chain is the longest. Honest blocks will rarely be part of the main chain.
    if beta >= 0.5:
        chain_quality = 0.0
        print("For β >= 0.5, the adversary's chain will grow faster indefinitely.")
        print("The expected chain quality, therefore, approaches 0.")
    else:
        # Based on a Markov chain model of the selfish mining process,
        # the expected chain quality can be calculated with the following formula:
        # Chain Quality = [(1-β)(1-2β)(1+2β-pβ)] / [1-β-β²]

        numerator = (1 - beta) * (1 - 2 * beta) * (1 + 2 * beta - p * beta)
        denominator = 1 - beta - beta**2
        
        if denominator == 0:
             # This condition is not met for β < 0.5, so it's a theoretical safeguard.
            print("Error: Calculation resulted in division by zero.")
            return

        chain_quality = numerator / denominator

    # Output the detailed calculation steps as requested
    print("The formula for expected chain quality (CQ) for β < 0.5 is:")
    print("CQ = [(1 - β)(1 - 2β)(1 + 2β - pβ)] / [1 - β - β²]\n")
    print("Plugging in the given values:")
    
    # Numerator calculation
    print("Numerator calculation:")
    print(f"Numerator = (1 - {beta}) * (1 - 2*{beta}) * (1 + 2*{beta} - {p}*{beta})")
    print(f"          = {1-beta:.3f} * {1-2*beta:.3f} * {1 + 2*beta - p*beta:.3f}")
    print(f"          = {numerator:.4f}\n")

    # Denominator calculation
    print("Denominator calculation:")
    print(f"Denominator = 1 - {beta} - {beta}²")
    print(f"            = 1 - {beta} - {beta**2:.2f}")
    print(f"            = {denominator:.4f}\n")

    # Final result
    print("Result:")
    print(f"Expected Chain Quality = {numerator:.4f} / {denominator:.4f}")
    print(f"                         = {chain_quality:.4f}")

    print("\n<<<" + str(chain_quality) + ">>>")

# Execute the calculation
calculate_chain_quality()