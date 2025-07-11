import sys

def find_animal_that_hears_muscle_twitches():
    """
    Analyzes which animal from a list can hear the low-frequency
    sound of a human muscle twitch.
    """
    # Step 1: Define the frequency of a human muscle twitch in Hertz.
    # This is a low-frequency sound. We'll use 25 Hz as a representative value.
    muscle_twitch_freq_hz = 25
    print(f"Analyzing the ability to hear a muscle twitch at {muscle_twitch_freq_hz} Hz.\n")

    # Step 2: Define the approximate hearing range (min_hz, max_hz) for each animal.
    animal_hearing_ranges = {
        'A. Dog': (67, 45000),
        'B. Bat': (20000, 200000),
        'C. Herring': (1, 4000),      # Herrings are exceptionally sensitive to low frequencies (infrasound).
        'D. Whale': (15, 100000),     # Baleen whales, in particular, hear very low frequencies.
        'E. Human': (20, 20000)       # The twitch is at the very edge of or below human hearing.
    }

    print("Checking which animal's hearing range includes the target frequency...")
    potential_hearers = []
    for animal, (min_freq, max_freq) in animal_hearing_ranges.items():
        # Step 3: Check if the muscle twitch frequency is within the animal's range.
        if min_freq <= muscle_twitch_freq_hz <= max_freq:
            print(f"- {animal}: YES (Range: {min_freq} Hz to {max_freq} Hz)")
            potential_hearers.append(animal)
        else:
            print(f"- {animal}: NO (Range: {min_freq} Hz to {max_freq} Hz)")

    # Step 4: Conclude based on the analysis.
    print("\n--- Conclusion ---")
    print(f"Based on the data, several animals' hearing ranges technically include {muscle_twitch_freq_hz} Hz.")
    print("However, the Herring is famously known for its ability to detect very low-frequency sounds.")
    print("It has a specialized connection between its swim bladder and inner ear that makes it extremely sensitive to subtle pressure changes and sounds in the water, such as those made by a predator's moving muscles. Therefore, it is the best answer.")


find_animal_that_hears_muscle_twitches()