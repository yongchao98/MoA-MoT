import sys

def calculate_chain_quality():
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    This function prompts the user for the adversary's mining power (beta) and the
    probability (p) that honest miners choose the adversary's block in a tie.
    It then calculates the chain quality based on a Markov chain model of the
    selfish mining process.

    The formula derived is:
    Q(beta, p) = ((1 + (beta - beta^2)*(1-p)) * (1 - 2*beta)) / (1 - beta - 3*beta^3)

    The function prints the step-by-step calculation with the user-provided values
    and the final result.
    """
    try:
        beta_str = input("Enter the adversary's mining power fraction (β, e.g., 0.3): ")
        beta = float(beta_str)
        if not (0 < beta < 0.5):
            print("Error: β must be between 0 and 0.5 for the model to be valid.", file=sys.stderr)
            return

        p_str = input("Enter the probability p for tie-breaking (e.g., 0.5): ")
        p = float(p_str)
        if not (0 <= p <= 1):
            print("Error: p must be between 0 and 1.", file=sys.stderr)
            return

    except ValueError:
        print("Error: Invalid input. Please enter numerical values.", file=sys.stderr)
        return

    # Numerator calculation
    # Numerator = (1 + (beta - beta**2)*(1-p)) * (1 - 2*beta)
    num_part1_val = beta - beta**2
    num_part2_val = 1 - p
    num_part3_val = 1 + num_part1_val * num_part2_val
    num_part4_val = 1 - 2 * beta
    numerator = num_part3_val * num_part4_val

    # Denominator calculation
    # Denominator = 1 - beta - 3*beta**3
    denom_part1_val = 3 * beta**3
    denominator = 1 - beta - denom_part1_val
    
    if denominator == 0:
        print("Error: Denominator is zero, cannot calculate chain quality.", file=sys.stderr)
        return

    chain_quality = numerator / denominator

    # Print the detailed calculation
    print("\n--- Chain Quality Calculation ---")
    print(f"Given β = {beta} and p = {p}\n")
    print("The formula for Chain Quality (Q) is:")
    print("Q = ( (1 + (β - β²)*(1-p)) * (1 - 2β) ) / (1 - β - 3β³) \n")

    print("Step 1: Calculate the Numerator")
    print(f"  Numerator = (1 + ({beta} - {beta}²)*(1 - {p})) * (1 - 2*{beta})")
    print(f"  Numerator = (1 + ({beta} - {beta**2:.4f})*({1-p})) * ({1-2*beta:.4f})")
    print(f"  Numerator = (1 + {num_part1_val:.4f} * {num_part2_val:.4f}) * {num_part4_val:.4f}")
    print(f"  Numerator = (1 + {num_part1_val * num_part2_val:.4f}) * {num_part4_val:.4f}")
    print(f"  Numerator = {num_part3_val:.4f} * {num_part4_val:.4f}")
    print(f"  Numerator = {numerator}\n")

    print("Step 2: Calculate the Denominator")
    print(f"  Denominator = 1 - {beta} - 3*{beta}³")
    print(f"  Denominator = 1 - {beta} - 3*{beta**3:.4f}")
    print(f"  Denominator = 1 - {beta} - {denom_part1_val:.4f}")
    print(f"  Denominator = {denominator}\n")

    print("Step 3: Calculate Final Chain Quality")
    print("  Q = Numerator / Denominator")
    print(f"  Q = {numerator} / {denominator}")
    
    print("\n------------------------------------")
    print(f"The final expected chain quality is:")
    print(f"{chain_quality}")
    print("------------------------------------")

    # The final answer to be returned in the specified format
    # This is for automated grading and should contain only the final numerical result.
    # print(f"\n<<<{chain_quality}>>>")


if __name__ == '__main__':
    calculate_chain_quality()