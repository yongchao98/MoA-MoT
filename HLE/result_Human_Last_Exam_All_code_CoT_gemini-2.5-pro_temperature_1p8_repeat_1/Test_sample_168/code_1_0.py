def solve_animal_hearing_puzzle():
    """
    Analyzes animal hearing ranges to determine which might hear a human muscle twitch.
    """
    # Step 1: Define the frequency of a human muscle twitch.
    # Muscle contractions produce very low-frequency sounds (mechanomyography).
    muscle_twitch_frequency_min = 5  # Hz
    muscle_twitch_frequency_max = 25 # Hz
    print(f"The sound from a human muscle twitch is very low frequency.")
    print(f"Frequency Range of Muscle Twitch: {muscle_twitch_frequency_min} Hz to {muscle_twitch_frequency_max} Hz.\n")

    # Step 2: Compare this frequency with the hearing ranges of the animals.
    print("Analyzing the hearing ranges of the potential animals:")
    print(" - A. Dog: ~67 Hz to 45,000 Hz. Cannot hear below 67 Hz.")
    print(" - B. Bat: ~2,000 Hz to 110,000 Hz (ultrasonic). Cannot hear low frequencies.")
    print(" - C. Herring: Has a very wide range, from <1 Hz to over 3,000 Hz. Can hear infrasound.")
    print(" - D. Whale (Baleen): ~10 Hz to 39,000 Hz. Can hear some, but Herring are notably sensitive.")
    print(" - E. Human: ~20 Hz to 20,000 Hz. Can barely hear the upper range of a twitch, missing most of it.\n")

    # Step 3: Conclude which animal's range fits.
    print("Conclusion:")
    print("The Herring's hearing range starts below 1 Hz, which fully covers the sound produced by a muscle twitch.")
    
    # Final 'Equation' / Comparison:
    herring_min_freq = 1
    print(f"\nFinal Comparison:")
    print(f"Muscle Twitch ({muscle_twitch_frequency_min}-{muscle_twitch_frequency_max} Hz) is within Herring's hearing range (<{herring_min_freq} Hz and up).")

solve_animal_hearing_puzzle()