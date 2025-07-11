import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability of choosing the adversary's block in a tie
                   when an honest miner finds the next block (0 <= p <= 1).
    """
    # Validate inputs
    if not (0 < beta < 1):
        print("Error: Adversary's mining power (β) must be between 0 and 1.", file=sys.stderr)
        return
    if not (0 <= p <= 1):
        print("Error: Probability (p) must be between 0 and 1.", file=sys.stderr)
        return

    # If beta is 0.5 or greater, the adversary will dominate the chain.
    if beta >= 0.5:
        print("For β >= 0.5, the adversary's chain will grow indefinitely.")
        print("The expected chain quality approaches 0.")
        print("\nFinal Equation:")
        print("Chain Quality = 0")
        print("\nFinal Answer:")
        print("<<<0>>>")
        return

    # Calculate the terms of the formula
    term1 = 1 - beta
    term2 = 1 - 2 * beta
    term3 = 1 + 2 * beta - p * beta
    numerator = term1 * term2 * term3
    denominator = 1 - beta - beta**2

    if denominator <= 0:
        # This case corresponds to when the honest chain can never catch up,
        # but for beta < 0.5, the denominator is always positive.
        # Golden ratio conjugate: (sqrt(5)-1)/2 approx 0.618
        # Denominator is zero when beta is (sqrt(5)-1)/2
        print("Error: The denominator is zero or negative, which shouldn't happen for β < 0.5.", file=sys.stderr)
        return

    chain_quality = numerator / denominator

    # Output the final equation with the numbers substituted
    print("Calculating Chain Quality for β = {} and p = {}".format(beta, p))
    print("\nFormula: Q = [(1-β) * (1 - 2*β) * (1 + 2*β - p*β)] / [1 - β - β^2]")
    print("\nFinal Equation:")
    print("Q = [({:.4f}) * ({:.4f}) * ({:.4f})] / [{:.4f}]".format(term1, term2, term3, denominator))
    print("Q = [{:.4f}] / [{:.4f}]".format(numerator, denominator))
    print("Q = {:.4f}".format(chain_quality))
    print("\nFinal Answer:")
    print("<<<{}>>>".format(round(chain_quality, 4)))


if __name__ == '__main__':
    try:
        # Prompt user for beta and p
        beta_str = input("Enter the adversary's mining power portion (β), e.g., 0.3: ")
        beta_val = float(beta_str)
        p_str = input("Enter the tie-breaking probability (p), e.g., 0.5: ")
        p_val = float(p_str)
        print("-" * 20)
        calculate_chain_quality(beta_val, p_val)
    except ValueError:
        print("Invalid input. Please enter numeric values for β and p.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
