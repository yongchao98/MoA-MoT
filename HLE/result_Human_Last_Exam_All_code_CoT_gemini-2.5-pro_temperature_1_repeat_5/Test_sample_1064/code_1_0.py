import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power controlled by the adversary (0 < beta < 1).
        p (float): The probability an honest miner chooses the adversary's block in a tie (0 <= p <= 1).

    Returns:
        float: The expected chain quality.
    """
    if not (0 < beta < 1):
        print("Error: Beta (β) must be between 0 and 1.", file=sys.stderr)
        return None
    if not (0 <= p <= 1):
        print("Error: Probability p must be between 0 and 1.", file=sys.stderr)
        return None

    alpha = 1 - beta

    # Calculate the numerator and denominator of the chain quality formula
    numerator = (alpha**2) * (1 + alpha * beta * (2 - p))
    denominator = 1 - beta**2 + beta**3
    
    if denominator == 0:
        print("Error: Calculation resulted in division by zero.", file=sys.stderr)
        return None

    chain_quality = numerator / denominator

    # Print the final equation with the calculated numbers
    print(f"Given β = {beta} and p = {p}:")
    print(f"The expected chain quality is calculated as Numerator / Denominator.")
    print(f"Numerator = (1-{beta})² * (1 + (1-{beta})*{beta}*(2-{p})) = {numerator}")
    print(f"Denominator = 1 - {beta}² + {beta}³ = {denominator}")
    print(f"Expected Chain Quality = {numerator} / {denominator} = {chain_quality}")
    
    return chain_quality

if __name__ == '__main__':
    # Example values for beta and p, as they were not provided in the prompt.
    # beta is the adversary's mining power. A value of 0.3 is a common example.
    beta_adversary = 0.3
    # p is the probability honest miners choose the adversary's block in a tie.
    # p=0.5 means they choose randomly.
    p_tie_break = 0.5
    
    calculate_chain_quality(beta_adversary, p_tie_break)
    # The final answer is the numerical value of the chain quality.
    # To conform to the required output format, we calculate it again and wrap it.
    alpha = 1 - beta_adversary
    numerator = (alpha**2) * (1 + alpha * beta_adversary * (2 - p_tie_break))
    denominator = 1 - beta_adversary**2 + beta_adversary**3
    final_answer = numerator / denominator
    print(f"\n<<<__{final_answer}__>>>")
