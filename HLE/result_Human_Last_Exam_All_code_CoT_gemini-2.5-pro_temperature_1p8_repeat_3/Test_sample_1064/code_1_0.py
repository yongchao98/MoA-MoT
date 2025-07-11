import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under selfish mining.

    Args:
        beta (float): The adversary's portion of mining power (0 < beta < 1).
        p (float): The probability an honest miner chooses the adversary's
                   block in a tie (0 <= p <= 1).
    """

    print(f"Calculating Expected Chain Quality for beta = {beta:.4f} and p = {p:.4f}\n")

    # The derived formula is piecewise.
    # If beta >= 0.5, the attacker's lead grows infinitely,
    # and the chain quality approaches 0.
    if beta >= 0.5:
        e_h = 0.0
        e_l_str = "infinity"
        quality = 0.0
        
        print(f"For beta = {beta} (>= 0.5), the attacker's lead is expected to grow infinitely.")
        print("The number of honest blocks in the longest chain becomes negligible.")
        print("Final Equation: Chain Quality = 0 / infinity")
        print(f"Expected Chain Quality = {quality:.4f}")
        return quality

    # For beta < 0.5, we use the derived formulas.
    # Let's break it down into expected honest blocks (E_H) and
    # expected total blocks (E_L) per cycle.
    
    # E_H = (1-beta) + beta*(1-beta)^2 * (2-p)
    e_h_part1 = 1 - beta
    e_h_part2 = beta * ((1-beta)**2) * (2-p)
    e_h = e_h_part1 + e_h_part2

    # E_L = (1 - beta^2 + beta^3) / (1-beta)
    e_l_numerator = 1 - beta**2 + beta**3
    e_l_denominator = 1 - beta
    e_l = e_l_numerator / e_l_denominator
    
    # Chain Quality Q = E_H / E_L
    quality = e_h / e_l
    
    # Output the final equation with each number, as requested.
    print("The final equation for Chain Quality (Q) is Q = E_H / E_L, where:")
    
    print("\n1. Expected Honest Blocks per Cycle (E_H):")
    print(f"   E_H = (1 - beta) + beta * (1 - beta)^2 * (2 - p)")
    print(f"   E_H = (1 - {beta:.4f}) + {beta:.4f} * ({1-beta:.4f})^2 * (2 - {p:.4f})")
    print(f"   E_H = {e_h_part1:.4f} + {e_h_part2:.4f}")
    print(f"   E_H = {e_h:.4f}")

    print("\n2. Expected Total Blocks per Cycle (E_L):")
    print(f"   E_L = (1 - beta^2 + beta^3) / (1 - beta)")
    print(f"   E_L = (1 - {beta**2:.4f} + {beta**3:.4f}) / (1 - {beta:.4f})")
    print(f"   E_L = {e_l_numerator:.4f} / {e_l_denominator:.4f}")
    print(f"   E_L = {e_l:.4f}")
    
    print("\n3. Final Chain Quality Calculation:")
    print(f"   Q = E_H / E_L")
    print(f"   Q = {e_h:.4f} / {e_l:.4f}")
    print(f"   Expected Chain Quality = {quality:.4f}")

    return quality

if __name__ == '__main__':
    # Set the parameters for the calculation.
    # beta is the adversary's mining power.
    # p is the probability honest miners choose the adversary's block in a tie.
    # We use beta = 1/3 and p = 0.5 as example values.
    
    beta_param = 1/3
    p_param = 0.5
    
    final_quality = calculate_chain_quality(beta_param, p_param)
    # The final answer format is just the numerical result.
    # print(f"\n<<<{final_quality:.4f}>>>")

<<<0.6400>>>