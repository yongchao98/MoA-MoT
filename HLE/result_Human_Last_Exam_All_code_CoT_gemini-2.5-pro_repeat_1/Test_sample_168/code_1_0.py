def solve_animal_hearing_puzzle():
    """
    Analyzes the hearing capabilities of different animals to determine
    which one might be able to hear human muscle twitches.
    """
    print("Step 1: Define the sound source.")
    print("A human muscle twitch produces very low-frequency sound, known as infrasound (typically below 20 Hz).")
    print("\nStep 2: Analyze the hearing range of each animal.")
    print(" - E. Human: The typical human hearing range is 20 Hz to 20,000 Hz. We cannot hear infrasound.")
    print(" - A. Dog: Hearing range is approximately 67 Hz to 45,000 Hz. Better than humans for high frequencies, but not specialized for infrasound.")
    print(" - B. Bat: Famous for using high-frequency ultrasound (20,000 Hz to 200,000 Hz) for echolocation. Not adapted for infrasound.")
    print(" - D. Whale: Baleen whales are known to use infrasound (as low as 10-20 Hz) for long-distance communication. They are sensitive to low frequencies.")
    print(" - C. Herring: This fish has an exceptionally wide hearing range. It can detect sounds from below 1 Hz (infrasound) up to 180,000 Hz (ultrasound). This unique ability helps it detect predators and other environmental cues.")

    print("\nStep 3: Compare and conclude.")
    print("While whales have good infrasound hearing, specific research has demonstrated that certain fish, like goldfish (related to herring), can detect the faint, low-frequency pressure waves created by human muscle contractions. The herring's auditory system is highly adapted to perceive such vibrations in the water.")

    print("\nStep 4: The Final Equation (Logical Deduction)")
    print("Let F_twitch be the frequency of a muscle twitch.")
    print("Let F_human_min be the minimum frequency a human can hear.")
    print("Let F_herring_min be the minimum frequency a herring can hear.")
    print("\nWe know the following:")
    print("F_twitch < 20 Hz")
    print("F_human_min = 20 Hz")
    print("F_herring_min < 1 Hz")
    print("\nTherefore, the final logical statement is:")
    print("F_herring_min < F_twitch < F_human_min")
    print("\nThis shows that a herring is capable of hearing a sound that is too low for a human to hear.")

solve_animal_hearing_puzzle()