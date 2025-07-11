def solve_chain_quality():
    """
    Calculates the expected chain quality for a selfish mining scenario.

    The scenario is defined by:
    - beta: The proportion of mining power held by the adversary.
    - p: The probability that an honest miner chooses the adversary's block in a tie.

    Chain Quality is defined as the long-term ratio of honest blocks in the longest
    chain to the total number of blocks in that chain.
    """

    # --- INPUTS ---
    # The portion of mining power of the adversary. This must be < 0.5.
    beta = 0.3
    # The probability of choosing the adversary's block in a tie.
    p = 0.5

    print(f"Calculating chain quality for beta = {beta} and p = {p}\n")

    # The formula for the expected chain quality (Q) is derived from a
    # Markov chain model of the selfish mining process.
    # Q = Numerator / Denominator

    # Calculate the numerator of the formula
    numerator = ((1 - beta)**2) * (1 + beta * (1 - beta) * (2 - p))

    # Calculate the denominator of the formula
    denominator = 1 - beta**2 + beta**3

    if denominator == 0:
        print("Error: Denominator is zero, cannot calculate chain quality.")
        return

    # Calculate the final chain quality
    chain_quality = numerator / denominator

    # --- OUTPUT ---
    # Print the numbers used in the final equation as requested.
    print(f"The equation for chain quality is Q = Numerator / Denominator")
    print(f"Numerator = ((1 - {beta})**2) * (1 + {beta} * (1 - {beta}) * (2 - {p})) = {numerator}")
    print(f"Denominator = 1 - {beta}**2 + {beta}**3 = {denominator}")
    print("-" * 30)
    print(f"The expected chain quality is: {chain_quality}")


solve_chain_quality()

<<<0.6947070707070707>>>