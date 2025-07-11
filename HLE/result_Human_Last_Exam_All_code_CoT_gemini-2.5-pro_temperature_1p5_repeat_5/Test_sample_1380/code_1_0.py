def find_bachelier_discovery():
    """
    This function explains the physical process Louis Bachelier connected to financial markets.
    """
    physical_process = "Brownian motion"
    explanation = (
        "In his 1900 dissertation 'Théorie de la Spéculation', Louis Bachelier "
        "modeled the random fluctuations of stock prices.\n"
        "He discovered that the mathematical principles governing these price movements "
        "were identical to those describing {process}.\n\n"
        "{process} is the random, erratic movement of microscopic particles suspended "
        "in a fluid (like pollen in water), resulting from collisions with the molecules of the fluid.\n\n"
        "Bachelier's work was the first to use this advanced stochastic model, now also known as a Wiener process, "
        "to describe financial markets, predating Einstein's famous work on the subject by five years."
    ).format(process=physical_process)

    print(explanation)

find_bachelier_discovery()