def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a selfish mining attack.

    Args:
        beta (float): The adversary's portion of the total mining power (0 < beta < 1).
        p (float): The probability that an honest miner chooses the adversary's block
                   in a tie (0 <= p <= 1).

    Returns:
        float: The expected chain quality.
    """
    if not (0 < beta < 1 and 0 <= p <= 1):
        raise ValueError("Inputs must be in the range 0 < beta < 1 and 0 <= p <= 1.")

    # If the adversary's power is 1/2 or more, their private chain grows faster
    # than the public chain indefinitely. The chain quality approaches 0.
    if beta >= 0.5:
        return 0.0

    # Calculate the expected number of honest blocks per mining round.
    # E[H] = (1-β) * 1 + β(1-β)² * (p*1 + (1-p)*2)
    # E[H] = (1-β) + β(1-β)² * (2-p)
    expected_honest_blocks = (1 - beta) + beta * (1 - beta)**2 * (2 - p)

    # Calculate the expected total blocks added to the chain per mining round.
    # This sums the length increases from all possible outcomes, weighted by probability.
    # E[L] = (1-β)*1 (honest win) + β(1-β)*2 (tie scenario) + Σ_{k=2 to ∞} [β^k(1-β) * k] (adv lead k>1)
    # The simplified closed-form expression is: 1 - β² + β / (1-β)
    expected_total_blocks = 1 - beta**2 + beta / (1 - beta)
    
    # Chain quality is the ratio of honest blocks to total blocks
    if expected_total_blocks == 0:
        return 0.0 # Should not be reachable for beta < 1

    chain_quality = expected_honest_blocks / expected_total_blocks
    
    return chain_quality

if __name__ == "__main__":
    # Example values for demonstration
    beta_adversary_power = 0.3
    p_tie_prob = 0.5
    
    # Perform the calculation
    quality = calculate_chain_quality(beta_adversary_power, p_tie_prob)
    
    # Output the result with the numbers used in the final equation
    print(f"Given an adversary's mining power (beta) of {beta_adversary_power} and a tie-breaking probability (p) of {p_tie_prob}:")
    print("-" * 80)
    
    # Print the components of the calculation
    numerator_calc = f"(1 - {beta_adversary_power}) + {beta_adversary_power} * (1 - {beta_adversary_power})**2 * (2 - {p_tie_prob})"
    denominator_calc = f"1 - {beta_adversary_power}**2 + {beta_adversary_power} / (1 - {beta_adversary_power})"
    
    num_val = (1 - beta_adversary_power) + beta_adversary_power * (1 - beta_adversary_power)**2 * (2 - p_tie_prob)
    den_val = 1 - beta_adversary_power**2 + beta_adversary_power / (1 - beta_adversary_power)

    print(f"Expected Honest Blocks per Round = {numerator_calc} = {num_val:.6f}")
    print(f"Expected Total Blocks per Round  = {denominator_calc} = {den_val:.6f}")
    print("-" * 80)
    print(f"Expected Chain Quality = {num_val:.6f} / {den_val:.6f}")
    print(f"\nFinal Expected Chain Quality: {quality:.6f}")
    
    # Return the final numerical answer in the required format
    print(f"\n<<<For beta={beta_adversary_power} and p={p_tie_prob}, the chain quality is {quality:.6f}>>>")
