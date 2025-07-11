import sys

def solve():
    """
    Calculates the expected chain quality in a Bitcoin system
    with a selfish miner.
    """
    if len(sys.argv) != 3:
        print("Usage: python your_script_name.py <beta> <p>")
        print("  <beta>: adversary's portion of mining power (e.g., 0.3)")
        print("  <p>: probability of choosing adversary's block in a tie (e.g., 0.5)")
        return

    try:
        beta = float(sys.argv[1])
        p = float(sys.argv[2])
    except ValueError:
        print("Error: beta and p must be numeric values.")
        return

    if not (0 < beta < 1):
        print("Error: Beta (β) must be between 0 and 1.")
        return
    if not (0 <= p <= 1):
        print("Error: Probability p must be between 0 and 1.")
        return

    # Calculate the rate of honest blocks being added to the chain (R_H)
    # R_H = (1-beta)^2 * (1 + beta*(1-beta)*(2-p))
    term1_rh = (1 - beta)**2
    term2_rh = (1 + beta * (1 - beta) * (2 - p))
    rate_honest = term1_rh * term2_rh

    # Calculate the rate of adversary blocks being added to the chain (R_A)
    # R_A = 2*beta^2*(1-beta)^2 + p*beta*(1-beta)^3 + 2*beta^2 - beta^3
    term1_ra = 2 * (beta**2) * ((1 - beta)**2)
    term2_ra = p * beta * ((1 - beta)**3)
    term3_ra = 2 * (beta**2)
    term4_ra = beta**3
    rate_adversary = term1_ra + term2_ra + term3_ra - term4_ra

    # Total rate of blocks added to the chain
    rate_total = rate_honest + rate_adversary

    # Chain quality is the ratio of honest blocks to total blocks
    if rate_total == 0:
        # This case should not happen for 0 < beta < 1
        chain_quality = 0
    else:
        chain_quality = rate_honest / rate_total

    print(f"Given adversary mining power (β) = {beta} and tie-breaking probability (p) = {p}:")
    print("-" * 30)
    print("Step 1: Calculate the rate of honest blocks added to the chain (R_H).")
    print(f"R_H = (1 - {beta})² * (1 + {beta} * (1 - {beta}) * (2 - {p}))")
    print(f"R_H = {term1_rh:.4f} * (1 + {beta*(1-beta):.4f} * {2-p})")
    print(f"R_H = {term1_rh:.4f} * {term2_rh:.4f}")
    print(f"R_H = {rate_honest:.6f}")
    print("-" * 30)

    print("Step 2: Calculate the rate of adversary blocks added to the chain (R_A).")
    print(f"R_A = 2*({beta})²*(1-{beta})² + {p}*{beta}*(1-{beta})³ + 2*({beta})² - ({beta})³")
    print(f"R_A = {term1_ra:.4f} + {term2_ra:.4f} + {term3_ra:.4f} - {term4_ra:.4f}")
    print(f"R_A = {rate_adversary:.6f}")
    print("-" * 30)

    print("Step 3: Calculate the total chain quality.")
    print("Chain Quality = R_H / (R_H + R_A)")
    print(f"Chain Quality = {rate_honest:.6f} / ({rate_honest:.6f} + {rate_adversary:.6f})")
    print(f"Chain Quality = {rate_honest:.6f} / {rate_total:.6f}")
    print(f"Expected Chain Quality = {chain_quality:.6f}")
    print("-" * 30)
    
    # The final answer format as requested.
    # To run, save the code as a python file (e.g., quality.py) and execute from the terminal:
    # python quality.py <beta> <p>
    # e.g., python quality.py 0.3 0.5
    final_answer = f"<<<{chain_quality:.6f}>>>"

if __name__ == '__main__':
    solve()