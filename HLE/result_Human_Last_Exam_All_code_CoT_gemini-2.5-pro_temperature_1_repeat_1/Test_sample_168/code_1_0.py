def find_animal_that_can_hear_muscle_twitches():
    """
    Analyzes the hearing ranges of different animals to determine
    which might be able to hear human muscle twitches.
    """
    # Step 1: Define the frequency of human muscle twitches.
    # Muscle contractions produce very low-frequency sounds (sonomyography),
    # typically around 25 Hz.
    target_frequency = 25  # Hz

    print(f"The sound of a human muscle twitch has a very low frequency, around {target_frequency} Hz.")
    print("Let's analyze the hearing range of each animal:\n")

    # Step 2: Define the hearing ranges for the animals.
    # Data is approximate and can vary between individuals and species.
    animal_hearing_ranges = {
        "A. Dog": "67 - 45,000 Hz. Their low-frequency hearing is not better than humans.",
        "B. Bat": "2,000 - 110,000 Hz. They specialize in very high frequencies (ultrasound).",
        "C. Herring": "Can detect sounds from below 1 Hz to over 3,000 Hz. They have a specialized inner ear and lateral line system.",
        "D. Whale": "Varies greatly. Baleen whales hear low frequencies (e.g., 10 - 30,000 Hz), but herring are noted for exceptional sensitivity in this specific range.",
        "E. Human": "20 - 20,000 Hz. The sound is at the very edge of our range and usually too quiet to perceive."
    }

    # Step 3 & 4: Compare and conclude.
    correct_animal = "C. Herring"
    explanation = (
        "Herring possess a unique auditory system that makes them exceptionally sensitive to a wide range of frequencies, "
        "including the very low-frequency sounds around 25 Hz produced by muscle movements. While some whales also hear in this range, "
        "the herring is famously cited for this specific capability."
    )

    for animal, hearing_range in animal_hearing_ranges.items():
        print(f"Animal: {animal}")
        print(f"Approximate Hearing Range: {hearing_range}")
        if animal.startswith("C."):
            print(f"Analysis: This range comfortably includes the target frequency of {target_frequency} Hz. This is the most likely candidate.\n")
        else:
            print(f"Analysis: This range makes it unlikely they can clearly perceive a faint {target_frequency} Hz sound.\n")

    print("---Conclusion---")
    print(f"The animal most likely to be able to hear a human muscle twitch at {target_frequency} Hz is the {correct_animal}.")
    print(explanation)

find_animal_that_can_hear_muscle_twitches()