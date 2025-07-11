def find_common_modulation():
    """
    This function explains the reasoning behind identifying the most common modulation
    from Maqam Bayati on D from the given list.
    """
    maqam_bayati_d = {
        "tonic": "D",
        "second": "E (half-flat)",
        "third": "F",
        "ghammaz (fourth)": "G"
    }

    explanation = [
        "1. The base Maqam is Bayati on D. Its lower jins (scale fragment) consists of D, E half-flat, F, and G.",
        "2. A taqsim (improvisation) in this maqam often features modulations, which are temporary shifts to other jinsan or maqamat.",
        "3. We evaluate the given options based on common Arabic music performance practice.",
        "4. Most options (A, B, C, D, E, F, G) are highly unusual as they start on non-pivotal notes or use very distant jinsan.",
        "5. Options H and I are both common. Option I, Jins Saba on D, involves changing the home jins to the closely related Jins Saba.",
        "6. Option H, Jins Saba on E, refers to creating a Saba feel on the second degree of the scale (E half-flat).",
        "7. This specific modulation is a well-known, signature feature of the melodic development of Bayati, often called 'Saba Zamzama'.",
        "8. Because it's a classic and defining characteristic of a Bayati taqsim, it stands out as the most common and recognizable modulation from the list.",
    ]

    print("Step-by-step reasoning:")
    for step in explanation:
        print(f"- {step}")

    # The final answer is H
    final_answer = 'H'
    print(f"\nThe most common modulation listed is described in option H.")

find_common_modulation()
<<<H>>>