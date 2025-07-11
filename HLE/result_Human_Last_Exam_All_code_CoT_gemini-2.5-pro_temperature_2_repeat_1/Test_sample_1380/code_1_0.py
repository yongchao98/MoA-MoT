def solve_bachelier_riddle():
    """
    This function identifies and prints the answer to the historical question
    about Louis Bachelier's work.
    """
    # The physical process Bachelier identified as being mathematically identical
    # to the behavior of stock prices.
    physical_process = "Brownian motion"

    # A brief explanation of the connection.
    explanation = (
        f"The physical process that Louis Bachelier discovered was mathematically identical\n"
        f"to how stock option probabilities spread over time is known as {physical_process}.\n\n"
        f"In his 1900 thesis, Bachelier modeled the random, unpredictable fluctuations of the Paris stock market.\n"
        f"This was the first time that the mathematics for what we now call stochastic processes or 'random walks'—\n"
        f"which describe {physical_process} (the random movement of particles suspended in a fluid)—was applied to finance."
    )

    print(explanation)

solve_bachelier_riddle()