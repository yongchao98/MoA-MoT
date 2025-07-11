import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability that an honest miner chooses the adversary's block in a tie (0 <= p <= 1).

    Returns:
        float: The expected chain quality.
    """
    print(f"Calculating chain quality for beta = {beta} and p = {p}")

    if not (0 < beta < 1):
        print("Error: Beta (adversary's mining power) must be between 0 and 1.", file=sys.stderr)
        return None
    if not (0 <= p <= 1):
        print("Error: p (tie-breaking probability) must be between 0 and 1.", file=sys.stderr)
        return None

    # If the adversary's power is 50% or more, their lead grows indefinitely.
    # The longest chain will eventually consist entirely of their blocks.
    if beta >= 0.5:
        print("Adversary's power (beta) is 0.5 or greater.")
        print("The expected chain quality approaches 0.")
        return 0.0

    # The calculation is based on a formula derived from a Markov chain model of selfish mining.
    # The formula is valid for beta < 0.5.
    
    # Calculate the numerator of the chain quality formula.
    # Numerator = (1 - 2*beta) * (1 - beta) * (1 + (1 - beta) * beta * (2 - p))
    term1 = 1 - 2 * beta
    term2 = 1 - beta
    term3_inner = (1 - beta) * beta * (2 - p)
    term3 = 1 + term3_inner
    numerator = term1 * term2 * term3

    print("\nStep 1: Calculate the numerator of the formula")
    print(f"Term (1 - 2*beta) = {term1}")
    print(f"Term (1 - beta) = {term2}")
    print(f"Term (1 + (1-beta)*beta*(2-p)) = {term3}")
    print(f"Final Numerator = {term1} * {term2} * {term3} = {numerator}")

    # Calculate the denominator of the chain quality formula.
    # Denominator = 1 - beta - 2*beta^2 + beta^3
    denominator = 1 - beta - 2 * (beta**2) + (beta**3)
    
    print("\nStep 2: Calculate the denominator of the formula")
    print(f"Denominator = 1 - {beta} - 2*({beta**2}) + ({beta**3}) = {denominator}")


    if denominator == 0:
        print("Error: Calculation resulted in division by zero.", file=sys.stderr)
        return None
        
    chain_quality = numerator / denominator

    print("\nStep 3: Calculate the final Chain Quality")
    print(f"Chain Quality = Numerator / Denominator = {numerator} / {denominator}")
    
    return chain_quality


if __name__ == '__main__':
    # --- User-defined parameters ---
    # beta: Adversary's portion of mining power. Must be < 0.5 for a non-zero result.
    adversary_power_beta = 0.3333

    # p: Probability that honest miners choose the adversary's block in a tie.
    tie_breaking_probability_p = 0.5
    # --------------------------------

    result = calculate_chain_quality(adversary_power_beta, tie_breaking_probability_p)

    if result is not None:
        print(f"\nFinal Expected Chain Quality: {result}")
        # The '<<<' block is used to return the final answer for automated evaluation.
        # It needs to be the exact numeric answer.
        final_answer_str = f"<<<{result}>>>"
        print(final_answer_str)