import sys

def calculate_chain_quality():
    """
    Calculates the expected chain quality under a selfish mining attack.

    The user can modify the values for beta and p below.
    """
    # --- User-configurable parameters ---

    # beta (β): The adversary's portion of the total mining power.
    # Must be a value between 0 and 1.
    # Example: 0.33 represents 33% of the mining power.
    beta = 0.33

    # p: The probability that honest miners will choose the adversary's block
    # in case of a tie. Must be a value between 0 and 1.
    # Example: 0.5 means a 50% chance.
    p = 0.5

    # --- Input validation ---
    if not (0 <= beta < 1):
        print("Error: beta (β) must be in the range [0, 1).", file=sys.stderr)
        return
    if not (0 <= p <= 1):
        print("Error: p must be in the range [0, 1].", file=sys.stderr)
        return

    # --- Derivation and Calculation ---

    # The derived formula for expected chain quality (Q) is:
    # Q = ((1 - β)^2 * (1 + β * (1 - β) * (2 - p))) / (1 - β^2 + β^3)

    print("This script calculates the expected chain quality (Q) for given parameters.")
    print(f"Adversary mining power (β): {beta}")
    print(f"Honest miner tie-breaking probability (p): {p}\n")

    # To fulfill the requirement of showing each number, we will substitute
    # the values into the formula and print each step.

    print("Step 1: Write down the formula with the given values.")
    formula_str = f"Q = ((1 - {beta})^2 * (1 + {beta} * (1 - {beta}) * (2 - {p}))) / (1 - {beta}^2 + {beta}^3)"
    print(formula_str + "\n")

    # Perform intermediate calculations for clarity
    one_minus_beta = 1 - beta
    two_minus_p = 2 - p
    beta_sq = beta ** 2
    beta_cub = beta ** 3

    print("Step 2: Calculate the intermediate terms.")
    step2_str = f"Q = (({one_minus_beta:.4f})^2 * (1 + {beta} * {one_minus_beta:.4f} * {two_minus_p:.4f})) / (1 - {beta_sq:.4f} + {beta_cub:.4f})"
    print(step2_str + "\n")
    
    # Calculate the values of the larger components
    term_in_paren = beta * one_minus_beta * two_minus_p
    sq_term = one_minus_beta**2
    denominator = 1 - beta_sq + beta_cub

    print("Step 3: Simplify the expression further.")
    step3_str = f"Q = ({sq_term:.4f} * (1 + {term_in_paren:.4f})) / {denominator:.4f}"
    print(step3_str + "\n")

    # Calculate the final numerator
    numerator = sq_term * (1 + term_in_paren)
    
    print("Step 4: Calculate the final numerator and denominator.")
    step4_str = f"Q = {numerator:.6f} / {denominator:.6f}"
    print(step4_str + "\n")

    # Final calculation
    if denominator == 0:
        # This case is unlikely for valid beta but is good practice to handle.
        final_quality = float('inf')
    else:
        final_quality = numerator / denominator

    print("Step 5: The final result for the expected chain quality is:")
    print(f"Q = {final_quality}")
    
    # The final answer in the required format
    print(f"\n<<<{final_quality}>>>")

if __name__ == '__main__':
    calculate_chain_quality()
