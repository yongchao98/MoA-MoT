def get_bachelier_discovery_info():
    """
    This function provides information about Louis Bachelier's discovery.
    """
    process_name = "Brownian motion"
    explanation = (
        "The specific physical process Louis Bachelier discovered was mathematically identical\n"
        "to how stock option probabilities spread over time is Brownian motion.\n\n"
        "In his 1900 thesis, 'Théorie de la Spéculation,' Bachelier modeled the random\n"
        "fluctuations of the French stock market. He realized that this financial process\n"
        "followed a random walk, which is mathematically described by the same principles\n"
        "that govern the random, erratic movement of particles suspended in a fluid—a\n"
        "phenomenon known as Brownian motion."
    )
    print(explanation)

if __name__ == "__main__":
    get_bachelier_discovery_info()