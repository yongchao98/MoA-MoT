def find_animal_that_hears_muscle_twitch():
    """
    Analyzes the hearing ranges of different animals to determine which one
    might be able to hear a human muscle twitch.
    """
    # Approximate low-end of hearing frequency in Hertz (Hz) for each animal
    # Note: These are general estimates and can vary.
    animal_hearing_low_freq = {
        'A. Dog': 67,
        'B. Bat': 2000,
        'C. Herring': 1,  # Herring can hear infrasound
        'D. Whale': 10, # Baleen whales can hear very low frequencies
        'E. Human': 20
    }

    # The frequency range of a human muscle twitch is very low
    muscle_twitch_min_freq = 5
    muscle_twitch_max_freq = 50

    print(f"A human muscle twitch produces a very low-frequency sound, around {muscle_twitch_min_freq} Hz to {muscle_twitch_max_freq} Hz.")
    print("-" * 50)
    print("Comparing this to the low-frequency hearing limits of the animals:")

    best_candidate = None
    min_freq_heard = float('inf')

    for animal, low_freq in animal_hearing_low_freq.items():
        print(f" - {animal}: Can hear down to approximately {low_freq} Hz.")
        if low_freq <= muscle_twitch_max_freq:
            # Check for the animal that can hear the lowest frequencies within the relevant range.
            # The lower an animal can hear, the more likely it is to detect faint infrasound.
            if low_freq < min_freq_heard:
                min_freq_heard = low_freq
                best_candidate = animal

    print("-" * 50)
    print(f"Based on the data, the Herring is the best candidate. Its ability to hear frequencies as low as {min_freq_heard} Hz, combined with a highly sensitive auditory system, makes it the most likely animal on the list to be able to detect the faint, low-frequency sound of a muscle twitch.")


if __name__ == "__main__":
    find_animal_that_hears_muscle_twitch()