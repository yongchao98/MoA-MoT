def solve_animal_hearing_puzzle():
    """
    Analyzes and prints the reasoning for determining which animal can hear human muscle twitches.
    """
    # Step 1: Define the frequency of the sound source.
    # The 'equation' involves comparing the muscle twitch frequency to animal hearing ranges.
    muscle_twitch_min_freq = 5  # in Hertz
    muscle_twitch_max_freq = 35 # in Hertz

    print("The task is to find which animal can hear human muscle twitches.")
    print("This requires comparing the sound frequency of a muscle twitch to the hearing range of each animal.")
    print(f"The sound of a muscle twitch (phonomyography) is in the very low-frequency range of approximately {muscle_twitch_min_freq} Hz to {muscle_twitch_max_freq} Hz.")
    print("-" * 30)

    # Step 2: Define the hearing ranges of the animals.
    # This data includes all the numbers needed for our comparison.
    hearing_ranges = {
        "A. Dog": "67 - 45,000 Hz",
        "B. Bat": "2,000 - 110,000 Hz",
        "C. Herring": "<1 - 180,000 Hz",
        "D. Whale (Baleen)": "10 - 30,000 Hz",
        "E. Human": "20 - 20,000 Hz"
    }

    print("Approximate hearing ranges of the animals:")
    for animal, freq_range in hearing_ranges.items():
        print(f"{animal}: {freq_range}")
    print("-" * 30)

    # Step 3: Analyze and conclude.
    print("Analysis:")
    print(f"The muscle twitch range ({muscle_twitch_min_freq}-{muscle_twitch_max_freq} Hz) is below or at the lowest limit for Humans and Dogs, and far too low for Bats.")
    print("While a Whale's hearing extends into this low-frequency range, the Herring is famously known for its extraordinarily wide and sensitive hearing.")
    print("A Herring's specialized auditory system, connecting its swim bladder to its inner ear, makes it extremely sensitive to pressure waves in water, allowing it to detect such faint, low-frequency sounds.")
    print("\nConclusion: The Herring is the animal most likely able to hear human muscle twitches.")

solve_animal_hearing_puzzle()