def trace_word_form():
    """
    Traces a hypothetical PIE root to a Middle English verb form
    and explains each step of the derivation.
    """
    pie_root = "*kʷeys-"
    pie_meaning = "to see, to heed"
    
    print(f"This script will trace the PIE root {pie_root} ({pie_meaning}) to a hypothetical Middle English verb form.\n")

    # Step 1: Proto-Indo-European to Proto-Germanic
    print("--- Step 1: PIE to Proto-Germanic (PGmc) ---")
    print(f"1. Start with the PIE root: {pie_root}")
    print("2. Form an o-grade causative. The e-grade '*e' becomes o-grade '*o'. The stem becomes *kʷoys-.")
    print("3. Add the causative suffix, which results in a Class 1 weak verb in Proto-Germanic with the infinitive ending *-janą.")
    print("4. Apply Grimm's Law: PIE *kʷ -> PGmc *hʷ.")
    print("5. Apply vowel change: PIE diphthong *oy -> PGmc diphthong *ai.")
    pgmc_form = "*hʷaisijaną"
    print(f"Resulting Proto-Germanic infinitive: {pgmc_form} (approx. 'to show, to inform')\n")
    
    # Step 2: Proto-Germanic to Old English (OE)
    print("--- Step 2: PGmc to Old English (OE) ---")
    print(f"1. The PGmc verb {pgmc_form} is a Class 1 weak verb with a long stem.")
    print("2. The '-j-' in the suffix causes i-umlaut on the root vowel.")
    print("3. The vowel change is: PGmc *ai -> pre-OE *ā -> (i-umlaut) OE ǣ (/æː/).")
    print("4. The initial PGmc *hʷ becomes the written cluster 'hw'.")
    oe_infinitive = "hwǣsan"
    print(f"Resulting Old English infinitive: {oe_infinitive}")
    print("5. To get the 3rd person singular present form ('he shows'), we add the ending '-þ' to the stem 'hwǣs-'.")
    oe_3sg_present = "hwǣsþ"
    print(f"Resulting Old English 3rd person singular present: {oe_3sg_present}\n")

    # Step 3: Old English to Middle English (ME)
    print("--- Step 3: OE to Middle English (ME) ---")
    print(f"1. Start with the OE form: {oe_3sg_present}")
    print("2. The initial OE cluster 'hw' is retained, with its spelling changed to 'wh'.")
    print("3. The OE long vowel 'ǣ' (/æː/) is raised to a long open-mid vowel 'ē' (/ɛː/), spelled 'e'.")
    print("4. The final OE consonant cluster '-sþ' simplifies to '-st' through assimilation.")
    me_form = "whest"
    print(f"Resulting Middle English form: {me_form}\n")
    
    # Final breakdown
    print("--- Final Form Breakdown ---")
    print(f"The final predicted Middle English form is: {me_form}")
    print("It is composed of the following parts based on sound changes:")
    print("wh- < OE hw- < PGmc *hʷ- < PIE *kʷ-")
    print("-e-  < OE ǣ < i-umlaut of *ā < PGmc *ai < PIE *oy")
    print("-st  < OE -sþ (from stem -s- + 3sg ending -þ)")

trace_word_form()
<<<whest>>>