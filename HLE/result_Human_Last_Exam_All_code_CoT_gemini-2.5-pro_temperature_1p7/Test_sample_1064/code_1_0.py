import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability that honest miners choose the adversary's block in a tie (0 <= p <= 1).
    """
    # Validate inputs
    if not (0 < beta < 1):
        print("Error: Adversary's mining power 'beta' must be between 0 and 1 (exclusive).", file=sys.stderr)
        return
    if not (0 <= p <= 1):
        print("Error: Probability 'p' must be between 0 and 1 (inclusive).", file=sys.stderr)
        return

    # If beta >= 0.5, the selfish miner's chain will always grow faster in expectation.
    # The attacker's lead becomes a transient random walk, meaning it will grow to infinity.
    # Thus, the proportion of honest blocks on the main chain will tend to 0.
    if beta >= 0.5:
        print("For beta >= 0.5, the adversary's chain is expected to grow indefinitely.")
        print("The long-term chain quality approaches 0.")
        print("\n<<<0.0>>>")
        return

    # Chain quality is defined as Rate_Honest / (Rate_Honest + Rate_Adversary)
    # We can calculate numbers N_H and N_A that are proportional to these rates.
    #
    # N_H is proportional to the rate of honest blocks being added to the chain.
    # N_A is proportional to the rate of adversary blocks being added to the chain.

    # Calculation for N_H
    # N_H = (1-β) * (1 + β(1-β)(2-p))
    nh_val = (1 - beta) * (1 + beta * (1 - beta) * (2 - p))

    # Calculation for N_A
    # N_A = (2β^2 - β^3)/(1-β) + β(1-β)(2β + p(1-β))
    na_term1 = (2 * beta**2 - beta**3) / (1 - beta)
    na_term2 = beta * (1 - beta) * (2 * beta + p * (1 - beta))
    na_val = na_term1 + na_term2

    # Total rate is proportional to N_H + N_A
    total_n = nh_val + na_val

    # Final chain quality
    if total_n == 0:
        # This case is unlikely with valid inputs but good practice to handle.
        chain_quality = 0
    else:
        chain_quality = nh_val / total_n
        
    # --- Outputting the result as requested ---
    
    print("Calculating the expected chain quality for a selfish miner.")
    print(f"Given parameters: beta = {beta}, p = {p}\n")

    print("The expected chain quality is calculated as: N_H / (N_H + N_A)\n")

    print("1. Calculating N_H (proportional to honest block rate):")
    print(f"   N_H = (1-{beta}) * (1 + {beta}*(1-{beta})*(2-{p}))")
    print(f"   N_H = {1-beta} * (1 + {beta*(1-beta)*(2-p):.4f})")
    print(f"   N_H = {nh_val:.4f}\n")

    print("2. Calculating N_A (proportional to adversary block rate):")
    print(f"   N_A = (2*{beta}² - {beta}³)/(1-{beta}) + {beta}*(1-{beta})*(2*{beta} + {p}*(1-{beta}))")
    print(f"   Term 1 = (2*{beta**2:.4f} - {beta**3:.4f})/{1-beta} = {na_term1:.4f}")
    print(f"   Term 2 = {beta*(1-beta):.4f} * (2*{beta} + {p*(1-beta):.4f}) = {na_term2:.4f}")
    print(f"   N_A = {na_term1:.4f} + {na_term2:.4f} = {na_val:.4f}\n")
    
    print("3. Final Chain Quality Calculation:")
    print(f"   Chain Quality = {nh_val:.4f} / ({nh_val:.4f} + {na_val:.4f})")
    print(f"   Chain Quality = {nh_val:.4f} / {total_n:.4f}")
    print(f"   Chain Quality = {chain_quality:.4f}\n")
    
    print(f"<<<{chain_quality:.6f}>>>")


if __name__ == '__main__':
    # You can change these values to test different scenarios
    # Adversary's mining power (e.g., 30%)
    adversary_power_beta = 0.3
    # Probability honest miners choose the adversary's block in a tie (e.g., 50%)
    tie_breaking_prob_p = 0.5
    
    calculate_chain_quality(adversary_power_beta, tie_breaking_prob_p)
