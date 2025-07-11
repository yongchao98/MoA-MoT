def find_bacheliers_process():
    """
    This function identifies and prints the name of the physical process
    Louis Bachelier used as a model for stock price fluctuations in his
    1900 thesis, "Théorie de la Spéculation".
    """
    # The physical process is the random movement of particles suspended in a fluid.
    # This was observed by botanist Robert Brown in 1827.
    # Bachelier was the first to create a mathematical model for this process
    # and apply it to financial markets, five years before Albert Einstein's famous paper on it.
    process_name = "Brownian motion"

    print("The specific physical process that Louis Bachelier discovered was mathematically identical to how stock option probabilities spread over time is:")
    print(process_name)

if __name__ == "__main__":
    find_bacheliers_process()