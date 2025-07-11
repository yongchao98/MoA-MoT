def analyze_moonlight_sonata():
    """
    This function provides a step-by-step music theory analysis
    for the provided excerpt of Beethoven's "Moonlight" Sonata.
    """

    # --- Part 1: Key of Modulation ---
    key_of_modulation = "F# minor"
    print("--- Question 1: To what key the music modulates in measure 11 and the beginning of measure 12? ---")
    print(f"Answer: The music modulates to (or more accurately, tonicizes) {key_of_modulation}.")
    print("Justification: The clearest evidence is on the third beat of measure 11, which features an E-sharp diminished triad (E#-G#-B). E-sharp is the leading tone to F-sharp, making this chord a leading-tone chord (vii°) that points strongly toward F# minor as a temporary tonic.\n")

    # --- Part 2: Connection and Justification ---
    home_key = "C# minor"
    relationship = "subdominant (iv)"
    print(f"--- Question 2: Given the 4 sharps key environment, what is the connection and / or justification for this particular modulation? ---")
    print(f"Answer: The home key of the sonata is {home_key}. The key of {key_of_modulation} is the {relationship} of {home_key}.")
    print("Justification: Tonicizing the subdominant is a very common compositional technique. It creates harmonic interest and movement away from the tonic area, often as a way to prepare for a larger modulation or a return to the tonic. Here, it is part of a transitional passage leading to the relative major key.\n")

    # --- Part 3: Roman Numeral Marking ---
    roman_numeral = "vii°/VI"
    chord_notes = "G#, B-natural, and D-natural"
    target_chord = "A major (the submediant, or VI chord, in C# minor)"
    print(f"--- Question 3: What would be the complete and accurate roman numeral marking for the first beat of measure 11? ---")
    print(f"Answer: The complete Roman numeral is {roman_numeral}.")
    print(f"Justification: The chord contains the notes {chord_notes}, forming a G-sharp diminished triad.")
    print(f"This chord functions as a secondary leading-tone chord. It is the leading-tone triad (vii°) of {target_chord}.")
    print("The numeral breakdown is as follows:")
    print("  vii°: A diminished triad built on the leading tone.")
    print("  /: Read as 'of', indicating it's a secondary function chord (i.e., not of the home tonic).")
    print(f"  VI: The chord being tonicized, which is the submediant of the home key, {home_key}.")

if __name__ == '__main__':
    analyze_moonlight_sonata()
