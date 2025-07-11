def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The adversary's portion of the mining power (0 < beta < 1).
        p (float): The probability of honest miners choosing the adversary's block in a tie (0 <= p <= 1).
    """
    print(f"Calculating for beta = {beta:.4f} and p = {p:.4f}:")
    print("-" * 40)

    if not (0 < beta < 1):
        print("Error: beta must be between 0 and 1.")
        return
    if not (0 <= p <= 1):
        print("Error: p must be between 0 and 1.")
        return

    # If beta >= 0.5, the selfish miner's chain will always grow faster in the long run.
    # The number of honest blocks in the longest chain will be negligible.
    if beta >= 0.5:
        chain_quality = 0
        print("Since beta >= 0.5, the selfish miner's chain will always be longer.")
        print(f"The expected chain quality approaches 0.")
        print(f"\nFinal Answer: {chain_quality:.4f}")
        return

    alpha = 1 - beta

    # Numerator of the formula: alpha * (alpha - beta) * (1 + 2*beta - beta*p)
    num_term1 = alpha
    num_term2 = alpha - beta
    num_term3 = 1 + 2 * beta - beta * p
    numerator = num_term1 * num_term2 * num_term3

    # Denominator of the formula: alpha^2 + alpha*beta - beta^2
    den_term1 = alpha**2
    den_term2 = alpha * beta
    den_term3 = beta**2
    denominator = den_term1 + den_term2 - den_term3

    if denominator == 0:
        print("Error: Calculation resulted in division by zero.")
        return
        
    chain_quality = numerator / denominator

    print(f"alpha = 1 - beta = 1 - {beta:.4f} = {alpha:.4f}\n")

    print("Numerator = alpha * (alpha - beta) * (1 + 2*beta - beta*p)")
    print(f"          = {num_term1:.4f} * ({num_term1:.4f} - {beta:.4f}) * (1 + 2*{beta:.4f} - {beta:.4f}*{p:.4f})")
    print(f"          = {num_term1:.4f} * {num_term2:.4f} * {num_term3:.4f}")
    print(f"          = {numerator:.4f}\n")

    print("Denominator = alpha^2 + alpha*beta - beta^2")
    print(f"            = {alpha:.4f}^2 + {alpha:.4f}*{beta:.4f} - {beta:.4f}^2")
    print(f"            = {den_term1:.4f} + {den_term2:.4f} - {den_term3:.4f}")
    print(f"            = {denominator:.4f}\n")

    print("Expected Chain Quality = Numerator / Denominator")
    print(f"                       = {numerator:.4f} / {denominator:.4f}")
    print(f"                       = {chain_quality:.4f}")
    
    # The final answer format requested by the user
    # print(f"\n<<< {chain_quality} >>>")


if __name__ == '__main__':
    # --- User-configurable parameters ---
    # Adversary's mining power (e.g., 0.25 for 25%, 1/3 for 33.3%)
    adversary_mining_power_beta = 1/3

    # Probability of honest miners choosing adversary's block in a tie
    tie_breaking_probability_p = 0.5
    # --- End of configuration ---

    calculate_chain_quality(adversary_mining_power_beta, tie_breaking_probability_p)
    
    # Example from the thought process to double check
    # print("\n" + "="*40 + "\n")
    # calculate_chain_quality(0.25, 0.5)