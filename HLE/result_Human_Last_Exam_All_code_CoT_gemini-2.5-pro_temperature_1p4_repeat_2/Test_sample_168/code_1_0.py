def find_animal_that_can_hear_muscle_twitches():
    """
    Analyzes which animal can hear the low-frequency sound of a human muscle twitch.
    """
    # Step 1: Store the approximate lowest hearing frequency (in Hz) for each animal.
    # Data is based on scientific studies of animal hearing ranges.
    animal_hearing_low_freq = {
        'A. Dog': 67,
        'B. Bat': 2000,
        'C. Herring': 1,  # Herring are known to be extremely sensitive to infrasound.
        'D. Whale': 15, # Varies, but some whales are excellent infrasound hearers.
        'E. Human': 20
    }

    # Step 2: Define the frequency of a human muscle twitch.
    # Muscle contractions produce infrasound, typically below 25 Hz.
    muscle_twitch_freq = 25

    print(f"The sound of a human muscle twitch is a very low frequency, under {muscle_twitch_freq} Hz.")
    print("Let's check which animal's hearing is sensitive enough to detect this.\n")

    winner = None

    # Step 3 & 4: Iterate and print the comparison for each animal.
    for animal, low_freq in animal_hearing_low_freq.items():
        # The 'equation' is the comparison of the animal's low frequency hearing limit
        # against the muscle twitch frequency.
        print(f"Checking {animal}:")
        print(f"The lowest frequency it can hear is {low_freq} Hz.")
        
        if low_freq < muscle_twitch_freq:
            print(f"Result: The equation '{low_freq} Hz < {muscle_twitch_freq} Hz' is True.")
            print(f"Therefore, {animal} might be able to hear a muscle twitch.\n")
            # We select the best candidate. Herring have famously sensitive low-frequency hearing.
            if winner is None or low_freq < animal_hearing_low_freq[winner]:
                 winner = animal
        else:
            print(f"Result: The equation '{low_freq} Hz < {muscle_twitch_freq} Hz' is False.")
            print(f"Therefore, {animal} likely cannot hear a muscle twitch.\n")
    
    print(f"Conclusion: Among the choices, the Herring ({winner}) has the most suitable hearing range, extending down to {animal_hearing_low_freq[winner]} Hz, making it the most likely to detect such a sound.")


find_animal_that_can_hear_muscle_twitches()