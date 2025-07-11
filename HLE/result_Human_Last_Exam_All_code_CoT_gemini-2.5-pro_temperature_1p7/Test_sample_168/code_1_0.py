def find_animal_with_sensitive_hearing():
    """
    Analyzes which animal can hear a very low-frequency sound like a muscle twitch.
    """
    # Step 1: Define the target sound frequency and animal hearing ranges.
    muscle_twitch_freq = "3-5 Hz"
    hearing_ranges = {
        "Human": "20 - 20,000 Hz",
        "Dog": "67 - 45,000 Hz",
        "Bat": "2,000 - 110,000+ Hz",
        "Whale (Baleen)": "7 - 35,000 Hz",
        "Herring (Shad)": "Below 1 Hz - 180,000 Hz"
    }

    print("The task is to find which animal might hear a human muscle twitch.")
    print(f"A human muscle twitch creates a very low frequency sound around {muscle_twitch_freq}.\n")

    print("Let's analyze the hearing ranges of the candidates:")
    for animal, freq_range in hearing_ranges.items():
        print(f"- {animal}: approx. {freq_range}")
    print("\n")

    # Step 2: Compare and reason.
    human_lower_limit = 20
    print(f"A Human's hearing does not go below {human_lower_limit} Hz, so we cannot hear our own muscles twitching.")
    print("A Dog's and Bat's hearing are specialized for higher frequencies, not infrasound.")
    print("While some Whales can hear infrasound down to 7 Hz, certain fish in the Herring family are known for extreme sensitivity to low-frequency pressure waves in water.")
    print("This sensitivity could allow them to detect the faint vibrations from a predator's muscle contractions.")
    print("\nTherefore, based on this biological data, the Herring is the animal famously cited for this ability.")

find_animal_with_sensitive_hearing()