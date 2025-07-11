import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a specific selfish mining strategy.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability of choosing the adversary's block in a tie
                     when an honest miner finds the tie-breaking block (0 <= p <= 1).

    Returns:
        float: The expected chain quality.
    """
    if not (0 < beta < 1):
        raise ValueError("Beta (β) must be between 0 and 1.")
    if not (0 <= p <= 1):
        raise ValueError("Probability p must be between 0 and 1.")

    alpha = 1 - beta

    # N_H is the expected number of honest blocks added to the chain per cycle.
    # A "cycle" is a sequence of events that starts in the normal state (state 0)
    # and ends when the system returns to it.
    n_h = alpha * (1 + beta * alpha * (2 - p))

    # n_total is the expected total number of blocks (honest + adversary) added
    # to the chain per cycle. This part of the formula is independent of p.
    n_total = (1 - beta**2 + beta**3) / alpha
    
    # n_a is the expected number of adversary blocks per cycle.
    n_a = n_total - n_h
    
    if n_total == 0:
        # Avoid division by zero, though this is only possible if n_h is also 0
        return 0.0
        
    chain_quality = n_h / n_total
    
    print("Inputs:")
    print(f"  Adversary mining power (β): {beta}")
    print(f"  Tie-breaking probability (p): {p}\n")
    
    print("Calculation:")
    print("The expected chain quality is the ratio of the expected number of honest blocks "
          "to the expected total number of blocks added to the chain per cycle.")
    print("\nExpected Honest Blocks per Cycle (N_H):")
    print(f"  N_H = (1 - {beta}) * (1 + {beta} * (1 - {beta}) * (2 - {p}))")
    print(f"  N_H = {alpha} * (1 + {beta * alpha} * {2-p})")
    print(f"  N_H = {n_h:.6f}\n")
    
    print("Expected Adversary Blocks per Cycle (N_A):")
    print(f"  N_A = {n_a:.6f}\n")

    print("Expected Total Blocks per Cycle (N_Total = N_H + N_A):")
    print(f"  N_Total = (1 - {beta}² + {beta}³) / (1 - {beta})")
    print(f"  N_Total = (1 - {beta**2} + {beta**3}) / {alpha}")
    print(f"  N_Total = {n_total:.6f}\n")

    print("Final Equation for Chain Quality (CQ):")
    print(f"  CQ = N_H / (N_A + N_H)")
    print(f"  CQ = {n_h:.6f} / ({n_a:.6f} + {n_h:.6f})")
    print(f"  CQ = {n_h:.6f} / {n_total:.6f}")
    print(f"  CQ = {chain_quality:.6f}")
    
    return chain_quality

if __name__ == '__main__':
    # --- You can change these values ---
    # beta is the adversary's fraction of the total mining power.
    beta_power = 0.35
    # p is the probability of adopting the adversary's block in a tie
    # when an honest miner finds the next block.
    tie_prob = 0.5
    # ---

    try:
        final_quality = calculate_chain_quality(beta_power, tie_prob)
        # The final answer is wrapped for the system to parse
        # This part will not be printed when you run the script
        # print(f"\n<<<ans:{final_quality:.6f}>>>")
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
