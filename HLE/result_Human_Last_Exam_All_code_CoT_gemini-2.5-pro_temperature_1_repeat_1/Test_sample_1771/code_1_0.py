def trace_etymology():
    """
    Traces the development of a PIE root through Germanic languages
    and prints the step-by-step derivation.
    """
    pie_verb = "*seh₂gieti"
    pgmc_verb = "*sōkijaną"
    on_verb = "sǿkja"
    oe_verb = "sēċan"
    ohg_verb = "suohhan"

    print(f"This analysis traces the PIE verb {pie_verb} to its descendants in several old Germanic languages.")
    print("The meaning evolved from 'giving a sign/tracking' to 'to seek'.\n")

    # --- Step 1: From Proto-Indo-European to Proto-Germanic ---
    print("1. Proto-Germanic Stage")
    print("-------------------------")
    print(f"The PIE root *seh₂g- develops into the Proto-Germanic (PGmc) verb {pgmc_verb} ('to seek').")
    print("Key sound changes:")
    print("  a) Grimm's Law: The PIE voiced stop *g becomes the PGmc voiceless stop *k.")
    print("  b) Vowel Merger: The o-grade of the root (*soh₂g-) becomes PGmc *sōk-, as the laryngeal *h₂ lengthens the vowel and then disappears.")
    print(f"Result in Proto-Germanic: {pgmc_verb}\n")

    # --- Step 2: From Proto-Germanic to Old Norse ---
    print("2. Old Norse Stage")
    print("--------------------")
    print(f"PGmc {pgmc_verb} becomes Old Norse (ON) {on_verb}.")
    print("Key sound change:")
    print("  a) i-Umlaut: The *j in the verb's suffix fronts the root vowel *ō to *ǿ.")
    print(f"Result in Old Norse: {on_verb}\n")

    # --- Step 3: From Proto-Germanic to Old English ---
    print("3. Old English Stage")
    print("----------------------")
    print(f"PGmc {pgmc_verb} becomes Old English (OE) {oe_verb}.")
    print("Key sound changes:")
    print("  a) i-Umlaut: The *j fronts the root vowel *ō to *ē.")
    print("  b) Palatalization: The consonant *k shifts to ċ (pronounced /tʃ/) before the new front vowel *ē.")
    print(f"Result in Old English: {oe_verb}\n")

    # --- Step 4: From Proto-Germanic to Old High German ---
    print("4. Old High German Stage")
    print("--------------------------")
    print(f"PGmc {pgmc_verb} becomes Old High German (OHG) {ohg_verb}.")
    print("Key sound changes:")
    print("  a) Diphthongization: The long vowel *ō breaks into the diphthong 'uo'.")
    print("  b) West Germanic Gemination: The consonant cluster *kj becomes a doubled fricative hh (pronounced /xx/).")
    print(f"Result in Old High German: {ohg_verb}\n")

    # --- Final Summary ---
    print("=== FINAL RESULTS ===")
    print(f"The PIE verb {pie_verb} yields the following words:")
    print(f"Proto-Germanic: {pgmc_verb}")
    print(f"Old Norse: {on_verb}")
    print(f"Old English: {oe_verb}")
    print(f"Old High German: {ohg_verb}")

if __name__ == '__main__':
    trace_etymology()