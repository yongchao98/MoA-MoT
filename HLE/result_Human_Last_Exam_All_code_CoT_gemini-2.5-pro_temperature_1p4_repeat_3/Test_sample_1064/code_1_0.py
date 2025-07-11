import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power controlled by the adversary (0 < beta < 1).
        p (float): The probability an honest miner chooses the adversary's block in a tie (0 <= p <= 1).
    
    Returns:
        float: The expected chain quality.
    """
    if not (0 < beta < 1):
        print("Error: beta (adversary's mining power) must be between 0 and 1.", file=sys.stderr)
        return None
    if not (0 <= p <= 1):
        print("Error: p (tie-breaking probability) must be between 0 and 1.", file=sys.stderr)
        return None

    # Expected number of honest blocks added to the chain per mining cycle.
    # This is (1-beta)*1 (for an honest block) + beta * (expected honest blocks from selfish attempt)
    honest_rate_numerator = 1 - beta
    honest_rate_selfish_part = beta * ((1 - beta)**2) * (2 - p)
    honest_rate = honest_rate_numerator + honest_rate_selfish_part

    # Expected total number of blocks added to the chain per mining cycle.
    # This is derived from analyzing the outcomes of all states.
    total_rate = 1 + beta + (beta**3) / (1 - beta)
    
    chain_quality = honest_rate / total_rate
    
    print("--- Calculating Expected Chain Quality ---")
    print(f"Given parameters: Adversary mining power (β) = {beta}, Tie-breaking probability (p) = {p}\n")
    
    print("Step 1: Calculate the effective rate of honest blocks being added to the main chain.")
    print("Formula: Rate(H) = (1-β) + β*(1-β)²*(2-p)")
    print(f"Calculation: Rate(H) = (1-{beta}) + {beta}*(1-{beta})²*(2-{p})")
    print(f"Calculation: Rate(H) = {honest_rate_numerator:.4f} + {honest_rate_selfish_part:.4f} = {honest_rate:.4f}\n")

    print("Step 2: Calculate the effective rate of total blocks being added to the main chain.")
    print("Formula: Rate(L) = 1+β + β³/(1-β)")
    print(f"Calculation: Rate(L) = 1+{beta} + {beta}³/(1-{beta})")
    print(f"Calculation: Rate(L) = {1+beta:.4f} + {beta**3:.4f}/{1-beta:.4f} = {total_rate:.4f}\n")
    
    print("Step 3: Calculate the Chain Quality.")
    print("Formula: Q = Rate(H) / Rate(L)")
    print(f"Calculation: Q = {honest_rate:.4f} / {total_rate:.4f}")
    print(f"Result: Expected Chain Quality = {chain_quality:.4f}\n")
    
    print("--- Final Algebraic Formula ---")
    print("Q = [ (1-β) + β(1-β)²(2-p) ] / [ 1+β + β³/(1-β) ]")

if __name__ == '__main__':
    # Example values for demonstration. You can change these.
    # β (beta) is the adversary's fraction of mining power.
    # p is the probability honest miners choose the selfish miner's block in a tie.
    beta_adversary = 0.3
    p_tie_break = 0.5
    
    calculate_chain_quality(beta_adversary, p_tie_break)
    
    # You can uncomment the line below to get the final formula as a string,
    # which is the formal answer to the question.
    final_formula = "[ (1-β) + β(1-β)²(2-p) ] / [ 1+β + β³/(1-β) ]"
    # print(f"\nFinal Answer: <<<{final_formula}>>>") # This is for generating the final output format.
    
# The final answer is the formula itself.
final_answer_formula = "[ (1-β) + β(1-β)²(2-p) ] / [ 1+β + β³/(1-β) ]"
print(f"\n<<<The expected chain quality is {final_answer_formula}>>>")
