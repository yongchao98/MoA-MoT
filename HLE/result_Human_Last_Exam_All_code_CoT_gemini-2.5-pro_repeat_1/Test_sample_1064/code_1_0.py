import argparse
import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The adversary's portion of mining power (0 <= beta <= 1).
        p (float): The probability an honest miner chooses the adversary's block in a tie (0 <= p <= 1).

    Returns:
        A tuple containing (Expected Honest Blocks, Expected Total Blocks, Chain Quality).
    """
    if not (0 <= beta <= 1):
        raise ValueError("Beta (β) must be between 0 and 1.")
    if not (0 <= p <= 1):
        raise ValueError("Probability p must be between 0 and 1.")

    # Handle edge cases to avoid division by zero and ensure correctness.
    if beta == 1.0:
        # If the adversary has all the mining power, no honest blocks will ever be in the longest chain.
        return 0.0, float('inf'), 0.0
    if beta == 0.0:
        # If the adversary has no mining power, all blocks in the longest chain will be honest.
        return 1.0, 1.0, 1.0

    alpha = 1.0 - beta

    # E_H is the expected number of honest blocks added to the main chain per renewal cycle.
    # The formula is derived from summing probabilities of different cycle paths.
    # Path 1 (Honest finds block): adds 1 honest block. Prob = alpha
    # Path 2 (Tie resolution): adds some honest blocks depending on p. Prob = beta*alpha
    expected_honest_blocks = alpha + beta * (alpha**2) * (2 - p)

    # E_L is the expected total number of blocks added to the main chain per renewal cycle.
    # This formula is also derived from the renewal cycle analysis.
    expected_total_blocks = (1 / alpha) - (beta**2)

    if expected_total_blocks == 0:
        # This case is unlikely with floating point arithmetic but good practice to handle.
        # If E_L is 0, and E_H is 0 (which it would be), quality is undefined.
        # However, for beta values approaching the root of 1/alpha - beta^2 = 0,
        # the chain quality approaches a specific limit. For simplicity, we can return 0.
        return expected_honest_blocks, expected_total_blocks, 0.0

    chain_quality = expected_honest_blocks / expected_total_blocks

    return expected_honest_blocks, expected_total_blocks, chain_quality

def main():
    """
    Main function to parse arguments and print the chain quality calculation.
    """
    parser = argparse.ArgumentParser(
        description="Calculate the expected chain quality under selfish mining.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('beta', type=float, help="Adversary's portion of mining power (e.g., 0.33)")
    parser.add_argument('p', type=float, help="Probability honest miners choose the adversary's block in a tie (e.g., 0.5)")

    if len(sys.argv) == 1:
        # Provide example usage if no arguments are given
        print("Usage: python your_script_name.py <beta> <p>")
        print("\nExample: Simulating a 1/3 power adversary (β=0.333) where network connectivity (p) is 0.5")
        print("python your_script_name.py 0.333 0.5")
        sys.exit(1)
        
    args = parser.parse_args()

    try:
        e_h, e_l, quality = calculate_chain_quality(args.beta, args.p)
        
        print(f"Given β = {args.beta} and p = {args.p}:\n")
        print("The expected chain quality is calculated as the ratio of expected honest blocks to expected total blocks added to the chain per cycle.\n")
        print(f"Expected Honest Blocks per Cycle (E_H) = {e_h}")
        print(f"Expected Total Blocks per Cycle (E_L)   = {e_l}")
        print("-" * 40)
        # The final equation with its components
        print(f"Expected Chain Quality = E_H / E_L")
        print(f"                       = {e_h} / {e_l}")
        print(f"                       = {quality}")

    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()