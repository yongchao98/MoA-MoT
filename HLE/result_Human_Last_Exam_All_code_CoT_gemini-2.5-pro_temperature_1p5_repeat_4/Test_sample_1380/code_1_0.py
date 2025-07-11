def explain_bachelier_discovery():
    """
    This function explains the physical process Louis Bachelier connected to financial markets.
    """
    process_name = "Brownian motion"
    explanation = f"""According to his 1900 PhD thesis, "Théorie de la Spéculation" (Theory of Speculation), Louis Bachelier discovered that the way stock and option price probabilities spread over time is mathematically identical to the physical process known as **{process_name}**.

{process_name} is the random movement of particles suspended in a fluid (a liquid or a gas). Bachelier's key insight was that the random, unpredictable fluctuations of the stock market could be modeled using the same stochastic mathematical framework that describes the diffusion of these particles. He was the first to apply this model to financial assets, predating Albert Einstein's famous 1905 paper on the same physical phenomenon.
"""
    print(explanation)

explain_bachelier_discovery()