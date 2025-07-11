import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power held by the adversary (0 <= beta < 1).
        p (float): The probability an honest miner chooses the adversary's block in a tie (0 <= p <= 1).
    """
    print(f"Calculating for beta = {beta} and p = {p}")

    # For beta >= 0.5, the selfish miner's lead grows indefinitely.
    # The probability of returning to a state where honest miners can contribute
    # becomes zero. The chain will be composed entirely of the adversary's blocks.
    if beta >= 0.5:
        print("Adversary's mining power (beta) is 0.5 or greater.")
        print("The adversary's chain will always be longer, leading to a chain quality of 0.")
        print("--- Final Equation ---")
        print("Number of Honest Blocks (proportional) = 0")
        print("Total Blocks in Longest Chain (proportional) = 1")
        print("\nExpected Chain Quality: 0.0")
        return

    # Based on the Markov chain model of selfish mining, the chain quality can be expressed as:
    # CQ = E[H] / E[L], where E[H] is the expected number of honest blocks and
    # E[L] is the expected total number of blocks added to the chain in one cycle.
    # After solving the state machine, we arrive at the following formula:
    
    # Numerator is proportional to the rate of honest blocks being added to the chain.
    numerator = ((1 - beta)**2) * (1 + beta * (1 - beta) * (2 - p))

    # Denominator is proportional to the total rate of growth of the chain.
    denominator = 1 - beta**2 + beta**3

    if denominator == 0:
        print("Error: Denominator is zero. Cannot calculate chain quality.")
        return

    chain_quality = numerator / denominator

    print("--- Final Equation ---")
    print(f"Number of Honest Blocks (proportional) = {numerator}")
    print(f"Total Blocks in Longest Chain (proportional) = {denominator}")
    print(f"Expected Chain Quality = {numerator} / {denominator}")
    print(f"\nFinal Result: {chain_quality}")


if __name__ == '__main__':
    # You can change these values to test different scenarios.
    # beta must be less than 0.5 for the selfish mining to not completely take over.
    try:
        # Check if command-line arguments are provided
        if len(sys.argv) == 3:
            b = float(sys.argv[1])
            prob_p = float(sys.argv[2])
        else:
            # Default values if no arguments are given
            print("Usage: python <script_name> <beta> <p>")
            print("Using default values: beta = 0.35, p = 0.5\n")
            b = 0.35
            prob_p = 0.5

        if not (0 <= b < 1):
            print("Error: beta must be between 0 and 1 (exclusive of 1).")
        elif not (0 <= prob_p <= 1):
            print("Error: p must be between 0 and 1.")
        else:
            calculate_chain_quality(b, prob_p)

    except ValueError:
        print("Error: Invalid input. Please provide numeric values for beta and p.")
    except Exception as e:
        print(f"An error occurred: {e}")

    # The expected chain quality for the derived formula
    final_beta = 0.35 if len(sys.argv) != 3 else float(sys.argv[1])
    final_p = 0.5 if len(sys.argv) != 3 else float(sys.argv[2])
    
    if 0 <= final_beta < 0.5 and 0 <= final_p <= 1:
        num = ((1 - final_beta)**2) * (1 + final_beta * (1 - final_beta) * (2 - final_p))
        den = 1 - final_beta**2 + final_beta**3
        result = num / den
        print(f"\n<<<Result for beta={final_beta}, p={final_p}: {result}>>>")
    elif final_beta >= 0.5 and 0 <= final_p <= 1:
        print(f"\n<<<Result for beta={final_beta}, p={final_p}: {0.0}>>>")
    else:
        print(f"\n<<<Result for beta={final_beta}, p={final_p}: Invalid input>>>")
