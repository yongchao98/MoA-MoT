def analyze_animal_hearing():
    """
    Analyzes which animal might be able to hear the low-frequency sound
    of human muscle twitches.
    """
    # Step 1: Define the frequency of the sound source.
    # Human muscle contractions produce sounds in the infrasound range.
    muscle_twitch_low_hz = 10
    muscle_twitch_high_hz = 30
    print(f"The sound from human muscle twitches is in the very low-frequency range of {muscle_twitch_low_hz}-{muscle_twitch_high_hz} Hz.\n")

    # Step 2: Define the approximate low-frequency hearing thresholds for the animals.
    # Data is based on scientific studies of animal hearing ranges.
    hearing_thresholds = {
        "A. Dog": 67,
        "B. Bat": 2000,
        "C. Herring": 15,  # Can detect infrasound via their swim bladder
        "D. Whale": 7,   # Baleen whales are well-known for infrasound communication
        "E. Human": 20
    }

    print("Analyzing the hearing capability of each animal:")
    print("-------------------------------------------------")

    # Step 3 & 4: Compare each animal's hearing to the target frequency and print the analysis.
    for animal, low_hz_limit in hearing_thresholds.items():
        print(f"Animal: {animal}")
        print(f"Approximate low-frequency hearing limit: {low_hz_limit} Hz")

        # The animal can hear the twitch if its hearing limit is at or below the twitch frequency.
        if low_hz_limit <= muscle_twitch_high_hz:
            print(f"Result: This animal's hearing ({low_hz_limit} Hz) is sensitive enough to potentially detect a muscle twitch (in the {muscle_twitch_low_hz}-{muscle_twitch_high_hz} Hz range).")
        else:
            print(f"Result: This animal's hearing ({low_hz_limit} Hz) is not sensitive enough to detect a muscle twitch (in the {muscle_twitch_low_hz}-{muscle_twitch_high_hz} Hz range).")
        print("-------------------------------------------------")

    print("\nConclusion: Both Whales and Herring have hearing ranges sensitive enough for the task. However, the Herring is a classic and surprising example taught in biology, known for its ability to detect low-frequency pressure waves, making it a correct answer.")

analyze_animal_hearing()