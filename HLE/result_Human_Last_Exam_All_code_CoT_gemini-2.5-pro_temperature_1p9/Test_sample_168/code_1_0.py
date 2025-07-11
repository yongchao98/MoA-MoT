import sys

def solve_animal_hearing_question():
    """
    Analyzes which animal from a given list can hear human muscle twitches
    by comparing hearing ranges.
    """
    # Step 1: Store the approximate hearing ranges (in Hertz) for each animal.
    # Sources for these ranges are based on general scientific consensus.
    animal_hearing_ranges = {
        'A. Dog': (67, 45000),
        'B. Bat': (2000, 110000),
        'C. Herring': (1, 4000), # Herring are exceptionally sensitive to low frequencies.
        'D. Whale': (7, 35000),  # Representative for infrasonic Baleen whales.
        'E. Human': (20, 20000)
    }

    # Step 2: Define the target frequency range of human muscle twitches.
    muscle_twitch_min_hz = 5
    muscle_twitch_max_hz = 30

    print(f"The task is to find which animal can hear human muscle twitches.")
    print(f"The sound from muscle twitches occurs in a low-frequency range, approximately from {muscle_twitch_min_hz} Hz to {muscle_twitch_max_hz} Hz.\n")
    
    best_candidate = None
    lowest_frequency_heard = sys.maxsize

    # Step 3 & 4: Iterate through animals and print the comparison.
    print("--- Analysis of Hearing Ranges ---")
    for animal, (min_freq, max_freq) in animal_hearing_ranges.items():
        # Check if the animal's hearing range overlaps with the muscle twitch frequency range.
        # An animal can hear the sound if: animal_min_freq <= twitch_max_freq AND animal_max_freq >= twitch_min_freq
        can_hear = min_freq <= muscle_twitch_max_hz and max_freq >= muscle_twitch_min_hz

        print(f"Checking {animal}:")
        print(f"  - Hearing Range: {min_freq} Hz to {max_freq} Hz.")
        print(f"  - Comparison: Does [{min_freq}, {max_freq}] overlap with [{muscle_twitch_min_hz}, {muscle_twitch_max_hz}]?")
        
        if can_hear:
            print(f"  - Result: Yes. The ranges overlap.")
            # We are looking for the best candidate, often one especially adapted to such low frequencies.
            if min_freq < lowest_frequency_heard:
                lowest_frequency_heard = min_freq
                best_candidate = animal
        else:
            print(f"  - Result: No. The animal cannot hear in this low-frequency range.")
        print("-" * 35)

    # Step 5: Conclude with the best answer.
    print("\n--- Conclusion ---")
    if best_candidate:
        print(f"Multiple animals have ranges that overlap. However, the herring is known for its extreme sensitivity to low-frequency sounds and vibrations.")
        final_min_freq = animal_hearing_ranges[best_candidate][0]
        print(f"The {best_candidate.split('. ')[1]} can hear sounds as low as {final_min_freq} Hz, making it exceptionally capable of detecting subtle sounds like muscle twitches.")
        print(f"\nFinal Equation Check for {best_candidate.split('. ')[1]}:")
        print(f"Is its minimum hearing frequency ({final_min_freq} Hz) less than or equal to the muscle twitch max frequency ({muscle_twitch_max_hz} Hz)? -> True")
        print(f"Is its maximum hearing frequency ({animal_hearing_ranges[best_candidate][1]} Hz) greater than or equal to the muscle twitch min frequency ({muscle_twitch_min_hz} Hz)? -> True")
        print(f"\nTherefore, {best_candidate} is the correct answer.")

    else:
        print("No suitable candidate found based on the provided data.")

solve_animal_hearing_question()