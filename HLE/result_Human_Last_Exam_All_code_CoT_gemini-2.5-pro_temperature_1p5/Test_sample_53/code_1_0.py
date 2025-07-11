import sys

def solve_historical_linguistics():
    """
    This script details the hypothetical derivation of a Middle English word
    from a Proto-Indo-European root, following standard sound changes.
    """

    # --- Stage 0: Define the forms at each historical stage ---
    pie_root = "*kʷeys-"
    pie_causative_stem = "*kʷoys-"
    proto_germanic_verb = "*hʷaizijaną"
    old_english_form = "hwǣrþ"
    middle_english_form = "whereth"

    print("This script calculates the hypothetical Middle English reflex of a PIE root.")
    print("-" * 70)

    # --- Stage 1: Proto-Indo-European to Proto-Germanic ---
    print("Step 1: From Proto-Indo-European to Proto-Germanic")
    print(f"   a. The starting PIE root is: {pie_root} ('to see, to heed')")
    print(f"   b. The o-grade causative stem is formed, changing the vowel to *o: {pie_causative_stem}")
    print("   c. Applying sound laws to create the Proto-Germanic verb:")
    print("      - Grimm's Law changes the initial *kʷ to *hʷ.")
    print("      - Verner's Law applies because the root syllable was unstressed, changing *s to *z.")
    print("      - The PIE diphthong *oy becomes Proto-Germanic *ai.")
    print("      - This creates the verb stem *hʷaiz-, which takes the Class 1 weak verb suffix *-janą.")
    print(f"   d. Resulting Proto-Germanic verb: {proto_germanic_verb}")
    print("-" * 70)

    # --- Stage 2: Proto-Germanic to Old English ---
    print("Step 2: From Proto-Germanic to Old English")
    print(f"   a. The starting Proto-Germanic form is: {proto_germanic_verb}")
    print("   b. Applying sound laws to create the Old English form:")
    print("      - In West Germanic, the *z between vowels becomes 'r' (Rhotacism).")
    print("      - The 'j' in the suffix causes i-umlaut, changing the stem diphthong *ai to Old English long *ǣ.")
    print("      - The 3rd person singular present ending for this verb class is -þ.")
    print(f"   c. Resulting Old English form ('he shows'): {old_english_form}")
    # Handle non-ASCII character for different terminals
    if sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
        print("      (Note: The form is 'hwaerth' using ASCII representation for 'ǣ')")
    print("-" * 70)

    # --- Stage 3: Old English to Middle English ---
    print("Step 3: From Old English to Middle English")
    print(f"   a. The starting Old English form is: {old_english_form}")
    print("   b. Applying sound laws to create the Middle English form:")
    print("      - The initial 'hw' sound comes to be spelled 'wh'.")
    print("      - The Old English long vowel *ǣ is raised and fronted, becoming the Middle English long open 'ē' sound, spelled 'e'.")
    print("      - The ending -þ evolves into the standard Middle English '-eth'.")
    print(f"   c. Resulting Middle English form: {middle_english_form}")
    print("-" * 70)

    # --- Final Summary ---
    print("Final Derivation Path (each 'number' in the 'equation'):")
    # This fulfills the request to "output each number in the final equation"
    final_equation = f"{pie_root} > {pie_causative_stem} > {proto_germanic_verb} > {old_english_form} > {middle_english_form}"
    print(final_equation)

solve_historical_linguistics()

<<<whereth>>>