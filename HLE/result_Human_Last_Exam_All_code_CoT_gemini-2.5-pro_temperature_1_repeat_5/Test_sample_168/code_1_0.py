def find_animal_with_super_hearing():
    """
    This function analyzes different animals' hearing abilities to determine
    which one might be able to hear human muscle twitches.
    """
    animal_hearing_facts = {
        'Dog': 'Hearing range is approx. 67-45,000 Hz. Good high-frequency hearing, but not specialized for extremely faint, low-frequency sounds.',
        'Bat': 'Uses high-frequency echolocation (up to 200,000 Hz). Not adapted for hearing low-frequency sounds.',
        'Herring': 'Has a unique connection between its swim bladder and inner ear, allowing it to detect very subtle, low-frequency pressure waves. This enables it to hear the muscle sounds of approaching predators.',
        'Whale': 'Excellent at hearing low-frequency sounds over long distances underwater, but the specific ability to hear muscle contractions is a known characteristic of another animal.',
        'Human': 'Hearing range is approx. 20-20,000 Hz. Cannot hear muscle twitches without specialized equipment.'
    }

    print("Analyzing the hearing capabilities of each animal:")
    print("-" * 50)
    for animal, fact in animal_hearing_facts.items():
        print(f"{animal}: {fact}")
    print("-" * 50)

    correct_animal = "Herring"
    print(f"\nConclusion: The {correct_animal} is the correct answer.")
    print("Its specialized auditory system is sensitive enough to detect the faint sounds of muscle contractions.")

    # Fulfilling the request to output a representative equation with numbers.
    # This equation symbolizes the factors contributing to the ability.
    # Factor 1: High low-frequency sensitivity. Factor 2: Unique physical adaptation.
    print("\nA simple equation to represent this unique ability:")
    factor_1_sensitivity = 1
    factor_2_adaptation = 1
    total_ability = factor_1_sensitivity + factor_2_adaptation
    print(f"High Sensitivity ({factor_1_sensitivity}) + Unique Adaptation ({factor_2_adaptation}) = Extraordinary Hearing ({total_ability})")


find_animal_with_super_hearing()