def bachelier_discovery():
    """
    This function explains the physical process Louis Bachelier linked to
    financial mathematics.
    """
    process_name = "Brownian motion"
    explanation = (
        "In his 1900 dissertation 'Théorie de la Spéculation', Louis Bachelier modeled the fluctuations\n"
        "of stock prices. He discovered that the mathematical model describing these price probabilities\n"
        "was identical to that of a well-known physical process: Brownian motion.\n\n"
        "Brownian motion is the random movement of particles (like dust or pollen) suspended in a fluid,\n"
        "caused by collisions with the fluid's molecules. Bachelier pioneered the mathematical theory for this\n"
        "process, which is now fundamental to both physics and financial modeling."
    )
    
    print(f"The physical process is: {process_name}")
    print("\n--- Explanation ---")
    print(explanation)

bachelier_discovery()