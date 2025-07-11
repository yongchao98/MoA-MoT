def find_animal_with_sensitive_hearing():
    """
    Analyzes the hearing capabilities of different animals to determine
    which might hear a human muscle twitch.
    """

    # A human muscle twitch produces a very faint sound around 25 Hz.
    muscle_twitch_hz = 25

    print(f"Analyzing which animal could hear a faint sound at ~{muscle_twitch_hz} Hz...")
    print("-" * 20)

    # Assigning a suitability score based on hearing capabilities.
    # Higher score means more likely to hear the sound.
    # 1: Unlikely, 5: Very Likely

    # Bat: Ultrasound specialist, does not hear low frequencies.
    bat_score = 1
    print("Bat: Hears very high frequencies (ultrasound). Score = 1")

    # Dog: Hears higher than humans, but not specialized for very low frequencies.
    dog_score = 2
    print("Dog: Hearing range starts around 67 Hz. Score = 2")

    # Human: The sound is at the absolute lower limit of our hearing.
    human_score = 3
    print("Human: Hearing range starts around 20 Hz, but sensitivity is very low. Score = 3")

    # Whale: Infrasound specialist, a good candidate.
    whale_score = 4
    print("Whale: Uses infrasound for communication. A strong candidate. Score = 4")

    # Herring: Has a specialized swim bladder-to-ear connection making it
    # extremely sensitive to faint, low-frequency pressure waves.
    herring_score = 5
    print("Herring: Unique anatomy makes it extremely sensitive to infrasound. The best candidate. Score = 5")
    print("-" * 20)

    # The user request requires outputting numbers from a final equation.
    # We will use the scores in a max() function as a representative equation.
    print("Final Equation (by suitability score):")
    print(f"max({bat_score}, {dog_score}, {human_score}, {whale_score}, {herring_score}) = {max(bat_score, dog_score, human_score, whale_score, herring_score)}")
    print("\nThe highest score corresponds to the Herring.")
    print("\nAnswer: C")

find_animal_with_sensitive_hearing()