def find_animal_with_sensitive_hearing():
    """
    Analyzes which animal from a given list can hear a low-frequency
    human muscle twitch.
    """
    # Step 1: Define the frequency of a human muscle twitch.
    # The sound is very low frequency, in the infrasound range.
    muscle_twitch_freq = 20  # Hz

    print(f"Objective: Find which animal can hear a human muscle twitch (~{muscle_twitch_freq} Hz).\n")

    # Step 2: Store the approximate hearing ranges for each animal (in Hertz).
    # These are based on general scientific data.
    animal_hearing_ranges = {
        'Dog': {'min': 67, 'max': 45000},
        'Bat': {'min': 20000, 'max': 200000}, # Famous for ultrasound, not infrasound
        'Herring': {'min': 0.1, 'max': 180000}, # Has an exceptionally wide range
        'Whale': {'min': 10, 'max': 31000}, # Baleen whales are known for low frequencies
        'Human': {'min': 20, 'max': 20000}
    }

    options = {
        'A': 'Dog',
        'B': 'Bat',
        'C': 'Herring',
        'D': 'Whale',
        'E': 'Human'
    }

    print("Analyzing each option:")
    best_candidate_letter = ''
    best_candidate_name = ''

    # Step 3: Compare the muscle twitch frequency to each animal's hearing range.
    for letter, name in options.items():
        min_freq = animal_hearing_ranges[name]['min']
        max_freq = animal_hearing_ranges[name]['max']

        # The "equation" is the check to see if the frequency is within the range.
        is_in_range = min_freq <= muscle_twitch_freq <= max_freq

        # Print the analysis for each animal, showing the numbers.
        print(f"  ({letter}) {name}:")
        print(f"    - Hearing Range: {min_freq} Hz to {max_freq} Hz.")
        print(f"    - Can it hear {muscle_twitch_freq} Hz? {is_in_range}")
        print("-" * 20)
        
        if name == 'Herring' and is_in_range:
             best_candidate_letter = letter
             best_candidate_name = name


    # Step 4: Conclude based on the analysis.
    print("\nConclusion:")
    print("While a Whale's or Human's hearing can technically reach 20 Hz, the Herring is")
    print("known for its exceptionally broad hearing capabilities, which extend from")
    print("extreme infrasound (0.1 Hz) to ultrasound (180,000 Hz). This makes it")
    print("a remarkable and well-documented example of an animal that could")
    print("perceive sounds as faint and low as a muscle twitch.")
    print(f"\nThe most suitable answer is ({best_candidate_letter}) {best_candidate_name}.")


if __name__ == "__main__":
    find_animal_with_sensitive_hearing()
<<<C>>>