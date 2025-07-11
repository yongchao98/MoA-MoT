import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a selfish mining attack.

    Args:
        beta (float): The adversary's portion of the total mining power (0 < beta < 0.5).
        p (float): The probability that honest miners choose the adversary's block in a tie (0 <= p <= 1).

    Returns:
        float: The expected chain quality.
    """
    if not (0 < beta < 1):
        print("Error: beta must be between 0 and 1.", file=sys.stderr)
        return None
    # For selfish mining to be profitable, beta should be > 1/3, and for the chain to not halt, beta < 0.5.
    # We will not enforce this here but it's a theoretical constraint.
    if beta >= 0.5:
        print("Warning: If beta >= 0.5, the adversary can trivially control the chain.", file=sys.stderr)
        
    if not (0 <= p <= 1):
        print("Error: p must be between 0 and 1.", file=sys.stderr)
        return None

    # Calculate the expected number of honest blocks added to the chain per cycle
    # E_H = (1-beta) + beta * [(1-beta)^2 * (2-p)]
    # E_H simplified:
    e_h = (1 - beta) + (2 - p) * beta * ((1 - beta) ** 2)

    # Calculate the expected total number of blocks added to the chain per cycle
    # E_T = (1-beta) + sum_{k=2 to inf}[beta^k*(1-beta)*k] + 2*beta^2*(1-beta) + 2*beta*(1-beta)^2
    # E_T simplified:
    e_t = (1 - beta**2 + beta**3) / (1 - beta)
    
    # Chain quality is the ratio of expected honest blocks to expected total blocks
    chain_quality = e_h / e_t
    
    print("Inputs:")
    print(f"  Adversary's mining power (β): {beta}")
    print(f"  Honest miners' preference in a tie (p): {p}")
    print("-" * 30)

    print("Calculation Steps:")
    # Show the calculation for E_H
    e_h_calc_str = f"E_H = 1 - {beta} + (2 - {p}) * {beta} * (1 - {beta})^2"
    print(e_h_calc_str)
    e_h_val_str = f"    = {1-beta} + {2-p} * {beta} * {round((1-beta)**2, 4)}"
    print(e_h_val_str)
    print(f"    = {e_h}")
    print()

    # Show the calculation for E_T
    e_t_calc_str = f"E_T = (1 - {beta}^2 + {beta}^3) / (1 - {beta})"
    print(e_t_calc_str)
    e_t_val_str = f"    = (1 - {round(beta**2, 4)} + {round(beta**3, 4)}) / {1-beta}"
    print(e_t_val_str)
    print(f"    = {e_t}")
    print()
    
    # Show the final calculation
    print("Final Chain Quality:")
    final_calc_str = f"Chain Quality = E_H / E_T = {e_h} / {e_t}"
    print(final_calc_str)
    print(f"              = {chain_quality}")
    
    return chain_quality

if __name__ == '__main__':
    # You can change these values to test different scenarios
    # beta: adversary's fraction of mining power.
    # A value of 1/3 is a critical threshold for selfish mining profitability.
    beta = 1/3
    
    # p: probability honest miners choose the attacker's block in a tie.
    # p=0.5 means they choose randomly.
    p = 0.5
    
    result = calculate_chain_quality(beta, p)
    # The final answer is wrapped in <<<>>>
    if result is not None:
        print(f"\n<<<{result}>>>")
