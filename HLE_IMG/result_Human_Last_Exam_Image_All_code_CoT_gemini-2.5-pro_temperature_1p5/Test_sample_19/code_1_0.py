import sys

def solve_music_theory():
    """
    This function analyzes the provided musical excerpt from Beethoven's Moonlight Sonata
    and answers the three-part question.
    """

    # --- Introduction and Plan ---
    print("Analyzing the musical excerpt to answer the three-part question.")
    print("The analysis will be broken down into three parts, followed by a final answer.\n")

    # --- Part 1: Key of Modulation ---
    print("--- Question 1: To what key does the music modulate in measures 11-12? ---\n")
    print("Step 1: The key signature has four sharps (F#, C#, G#, D#). This indicates either E Major or C-sharp minor.")
    print("Step 2: The piece clearly establishes the home key as C-sharp minor in the opening measures.")
    print("Step 3: In measures 10-12, the harmony points towards a new key. The chord on beat 1 of measure 11 is a B Major triad (notes: B, D#, F#).")
    print("Step 4: A B Major chord is the dominant (V) chord in the key of E Major.")
    print("Step 5: This V chord (B Major) strongly leads to the I chord (E Major) which arrives in measure 13, confirming the modulation.")
    print("\nAnswer 1: The music is modulating to the key of E Major.\n")

    # --- Part 2: Justification for Modulation ---
    print("--- Question 2: What is the connection/justification for this modulation? ---\n")
    print("Step 1: The original home key is C-sharp minor.")
    print("Step 2: The new key is E Major.")
    print("Step 3: E Major is the relative major of C-sharp minor. They are 'relative' because they share the exact same key signature (four sharps).")
    print("Step 4: In Western classical music, modulating from a minor key to its relative major is a very common and harmonically logical progression. It provides tonal contrast while maintaining a close relationship between the keys.")
    print("\nAnswer 2: The modulation is from the tonic minor (C-sharp minor) to its relative major (E Major), which is a standard and closely-related key change.\n")

    # --- Part 3: Roman Numeral Analysis ---
    print("--- Question 3: What is the Roman numeral for the first beat of measure 11? ---\n")
    print("Step 1: The chord is a B Major triad in root position.")
    print("Step 2: This chord functions as the dominant of the new key, E Major.")
    print("Step 3: For a complete analysis, we relate this function back to the original key of C-sharp minor.")
    print("Step 4: In C-sharp minor, the mediant (the third scale degree) corresponds to the 'III' chord, which is E Major.")
    print("Step 5: The B Major chord is the dominant of this III chord. This is called a 'secondary dominant'.")
    print("Step 6: This function is notated as 'V of III' (the dominant of the mediant).")
    print("\nAnswer 3: The complete and accurate Roman numeral marking is V/III.\n")

    # --- Final Combined Answer ---
    final_answer = (
        "1. The music modulates to E Major.\n"
        "2. The justification is that E Major is the relative major of the home key, C-sharp minor. This is a common and tonally close modulation because they share the same key signature.\n"
        "3. The Roman numeral marking is V/III. The chord is a B Major triad, which is the dominant (V) of E Major. Since E Major is the mediant (III) of the original key, the full analysis is 'V of III'."
    )
    
    # Python 3.6+ f-string is used here. For older versions, use .format()
    # or simple string concatenation.
    if sys.version_info >= (3, 6):
        sys.stdout.write(f"\n<<<3 parts answer:\n{final_answer}>>>\n")
    else:
        sys.stdout.write("\n<<<3 parts answer:\n" + final_answer + ">>>\n")

# Execute the analysis function
solve_music_theory()