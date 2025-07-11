import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The fraction of mining power of the adversary (0 <= beta < 1).
        p (float): The probability that honest miners choose the adversary's block in a tie (0 <= p <= 1).

    Returns:
        float: The expected chain quality, or None if inputs are invalid.
    """
    if not (0 <= beta < 1 and 0 <= p <= 1):
        print("Error: Invalid inputs. 'beta' must be in [0, 1) and 'p' in [0, 1].")
        return None

    # If beta >= 0.5, the selfish miner's chain grows indefinitely.
    # Honest blocks are rarely added to the longest chain. Chain quality approaches 0.
    if beta >= 0.5:
        print(f"For beta = {beta}, the selfish mining attack is fully effective.")
        print("Expected number of honest blocks on the main chain approaches 0.")
        print("Chain Quality = 0.0")
        return 0.0

    # Let H be the expected number of honest blocks and A be the expected number of
    # adversary blocks added to the main chain during one cycle starting from State 0.
    
    # Pre-calculate powers of beta
    beta_sq = beta ** 2
    beta_cb = beta ** 3
    one_minus_beta = 1 - beta
    one_minus_beta_sq = one_minus_beta ** 2

    # Expected number of honest blocks (H)
    # H = (1-beta) + beta * h_1
    # h_1 = (1-beta)^2 * (2-p)
    # H = (1-beta) + beta * (1-beta)^2 * (2-p)
    # After simplification:
    H = 1 + beta - 4 * beta_sq + 2 * beta_cb - p * beta * one_minus_beta_sq
    
    # Expected number of adversary blocks (A)
    # A = beta * a_1
    # a_1 = beta * (2 + beta/(1-beta)) + (1-beta) * (2*beta + p*(1-beta))
    # After simplification:
    A = 4 * beta_sq - 2 * beta_cb + beta_cb / one_minus_beta + p * beta * one_minus_beta_sq

    # The sum H+A simplifies to (1 - beta^2 + beta^3) / (1 - beta)
    total_blocks = (1 - beta_sq + beta_cb) / one_minus_beta
    
    # Chain quality is H / (H + A)
    if total_blocks == 0:
      # This case should not happen for valid beta < 0.5
      return 0.0

    chain_quality = H / total_blocks

    print("Calculation based on the standard selfish mining model:")
    print(f"Given beta = {beta} and p = {p}\n")
    print("Expected number of honest blocks in a cycle (H): {:.6f}".format(H))
    print("Expected number of adversary blocks in a cycle (A): {:.6f}".format(A))
    print("Expected total blocks in a cycle (H+A): {:.6f}".format(total_blocks))
    print("\nFinal Equation for Chain Quality:")
    print("Chain Quality = H / (H + A)")
    print("              = {:.6f} / {:.6f}".format(H, total_blocks))
    print("              = {:.6f}".format(chain_quality))
    
    return chain_quality

if __name__ == '__main__':
    # Example usage. You can change these values.
    # beta is the adversary's mining power portion.
    # p is the probability honest miners choose the adversary's block in a tie.
    try:
        # Check for command line arguments
        if len(sys.argv) == 3:
            beta_input = float(sys.argv[1])
            p_input = float(sys.argv[2])
        else:
            # Default values if no arguments are provided
            print("Usage: python your_script_name.py <beta> <p>")
            print("Using default values: beta = 0.3, p = 0.5\n")
            beta_input = 0.3
            p_input = 0.5

        result = calculate_chain_quality(beta_input, p_input)
        if result is not None:
            print(f"\n<<<{result:.6f}>>>")

    except (ValueError, IndexError):
        print("Error: Please provide valid numerical inputs for beta and p.")
        print("Usage: python your_script_name.py <beta> <p>")
