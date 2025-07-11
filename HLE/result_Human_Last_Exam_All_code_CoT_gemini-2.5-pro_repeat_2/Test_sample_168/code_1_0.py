def solve_animal_hearing_puzzle():
    """
    This script analyzes the hearing ranges of different animals to determine
    which might be able to hear the low-frequency sound of a human muscle twitch.
    """
    # Step 1: Define the target frequency range for a human muscle twitch.
    # This is our target "equation" to check against.
    twitch_min_hz = 5
    twitch_max_hz = 25

    # Step 2: Define the approximate hearing ranges for the given animals.
    hearing_ranges = {
        'Human': (20, 20000),
        'Dog': (67, 45000),
        'Bat': (20000, 200000),
        'Whale': (10, 30000), # Baleen whales are sensitive to very low frequencies.
        'Herring': (1, 3000)  # The clupeid family (herrings) can detect infrasound.
    }

    # Step 3: Print the analysis.
    print("Problem: Which animal might be able to hear human muscle twitches?")
    print("-" * 60)
    print("The frequency of a human muscle twitch needs to be compared against the hearing range of each animal.")
    # Here we output the numbers in the 'equation' as requested.
    print(f"Target Frequency Range of Muscle Twitch: {twitch_min_hz} Hz - {twitch_max_hz} Hz")
    print("-" * 60)
    
    best_candidate = None
    lowest_hearing_freq = float('inf')

    print("Analyzing each animal's hearing range:")
    for animal, (low_freq, high_freq) in hearing_ranges.items():
        # Output each number from the hearing range data.
        print(f"\nAnimal: {animal}")
        print(f"Hearing Range: {low_freq} Hz to {high_freq} Hz")
        
        # Check if the animal's hearing range overlaps with the muscle twitch frequency.
        if low_freq <= twitch_max_hz:
            print(f"Result: PLAUSIBLE. This animal can hear frequencies as low as {low_freq} Hz, which is in or below the muscle twitch range.")
            # Keep track of the animal that can hear the lowest frequencies.
            if low_freq < lowest_hearing_freq:
                lowest_hearing_freq = low_freq
                best_candidate = animal
        else:
            print(f"Result: UNLIKELY. The lower hearing limit of {low_freq} Hz is above the muscle twitch frequency range.")

    print("-" * 60)
    print("Conclusion:")
    print(f"The sound from muscle twitches is infrasound ({twitch_min_hz}-{twitch_max_hz} Hz).")
    print(f"The analysis shows that the '{best_candidate}' has the lowest hearing threshold ({lowest_hearing_freq} Hz), making it the most likely candidate to detect such a faint, low-frequency sound.")

# Execute the function to print the analysis.
solve_animal_hearing_puzzle()