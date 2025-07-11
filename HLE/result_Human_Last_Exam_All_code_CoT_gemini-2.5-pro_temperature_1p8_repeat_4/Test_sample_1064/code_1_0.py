def calculate_chain_quality():
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    This function prompts the user for the adversary's mining power (beta)
    and the tie-breaking probability (p), then computes and explains the
    resulting chain quality.
    """
    try:
        beta_str = input("Enter the adversary's mining power portion (β, a value between 0 and 1): ")
        beta = float(beta_str)
        if not (0 < beta < 1):
            print("Error: β must be strictly between 0 and 1.")
            return

        p_str = input("Enter the probability p for tie-breaking (a value between 0 and 1): ")
        p = float(p_str)
        if not (0 <= p <= 1):
            print("Error: p must be between 0 and 1 (inclusive).")
            return

    except ValueError:
        print("Invalid input. Please enter valid numbers.")
        return

    # --- Calculation ---
    # Expected number of honest blocks per cycle
    # E_H = (1-β) * [1 + β*(1-β)*(2-p)]
    e_h = (1 - beta) * (1 + beta * (1 - beta) * (2 - p))

    # Expected total number of blocks per cycle
    # E_L = 1 + β - β³
    e_l = 1 + beta - beta**3

    if e_l == 0:
        chain_quality = float('inf') if e_h != 0 else 0
    else:
        chain_quality = e_h / e_l

    # --- Output ---
    print("\n--------------------")
    print("The formula for Expected Chain Quality (CQ) is the ratio of expected honest blocks to expected total blocks per cycle:")
    print("CQ = E[Honest Blocks] / E[Total Blocks]")
    print("\nBased on the model, the formula is:")
    print("CQ = ( (1-β) * (1 + β*(1-β)*(2-p)) ) / ( 1 + β - β³ )")

    print(f"\nPlugging in your values β = {beta} and p = {p}:")

    # The components of the equation
    one_minus_beta = 1 - beta
    two_minus_p = 2 - p
    beta_cubed = beta**3
    
    # Print the equation with each number explicitly shown
    print(f"\nStep 1: Substitute β and p into the formula.")
    print(f"CQ = ( ({one_minus_beta:.4f}) * (1 + {beta:.4f} * {one_minus_beta:.4f} * {two_minus_p:.4f}) ) / ( 1 + {beta:.4f} - {beta_cubed:.4f} )")
    
    # Calculate intermediate values
    inner_term = beta * one_minus_beta * two_minus_p
    denominator = 1 + beta - beta_cubed
    
    print(f"\nStep 2: Calculate the products and sums inside the parentheses.")
    print(f"CQ = ( {one_minus_beta:.4f} * (1 + {inner_term:.4f}) ) / ( {denominator:.4f} )")

    # Calculate numerator
    numerator = one_minus_beta * (1 + inner_term)

    print(f"\nStep 3: Calculate the final numerator.")
    print(f"CQ = {numerator:.4f} / {denominator:.4f}")

    print("\n--------------------")
    print("Final Result:")
    print(f"The expected chain quality is: {chain_quality}")
    print("--------------------")

# Run the calculation
calculate_chain_quality()