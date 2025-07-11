def analyze_animal_hearing():
    """
    Analyzes which animal from a given list might be able to hear a human muscle twitch.
    """
    # Human muscle contractions produce a very low-frequency sound.
    muscle_twitch_freq_hz = 25

    # Approximate hearing ranges for the animals in Hz (min, max).
    # 'None' indicates no practical limit on that end for this context.
    animal_hearing_ranges = {
        "A. Dog": (67, 45000),
        "B. Bat": (20000, 200000),
        "C. Herring": (1, 3000), # Herring have a specialized wide range, very sensitive to low frequencies
        "D. Whale (Baleen)": (10, 30000),
        "E. Human": (20, 20000)
    }

    print(f"The task is to determine which animal can hear a sound at approximately {muscle_twitch_freq_hz} Hz.\n")

    for animal, (min_freq, max_freq) in animal_hearing_ranges.items():
        can_hear = "No"
        reason = ""
        if min_freq <= muscle_twitch_freq_hz <= max_freq:
            can_hear = "Yes (in theory)"
        
        if "Dog" in animal:
            reason = f"Their low-end hearing at ~{min_freq} Hz is too high."
        elif "Bat" in animal:
            reason = f"They specialize in ultrasound, hearing from {min_freq} Hz, far above the target."
        elif "Herring" in animal:
            reason = (f"Their hearing starts as low as {min_freq} Hz. They have specialized swim bladders "
                      "that detect low-frequency pressure waves, making them extremely sensitive to sounds "
                      "like a muscle twitch. This makes them the best candidate.")
        elif "Whale" in animal:
            reason = f"Their range starting at {min_freq} Hz includes the target frequency, making them a possible candidate."
        elif "Human" in animal:
            reason = (f"While {muscle_twitch_freq_hz} Hz is in our range (starting at {min_freq} Hz), the sound is too "
                      "faint for our ears to perceive audibly.")
        
        print(f"Animal: {animal}")
        print(f"Hearing Range: ~{min_freq} - {max_freq} Hz")
        print(f"Can hear {muscle_twitch_freq_hz} Hz? {can_hear}")
        print(f"Analysis: {reason}\n")
    
    print("Conclusion: Based on extreme low-frequency sensitivity, the Herring is the most likely animal to be able to hear a human muscle twitch.")

analyze_animal_hearing()