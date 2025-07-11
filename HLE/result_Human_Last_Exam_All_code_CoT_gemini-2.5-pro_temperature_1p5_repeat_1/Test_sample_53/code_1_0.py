def trace_word_evolution():
    """
    Traces a hypothetical PIE root to its Middle English form
    and prints each stage of the derivation.
    """
    pie_root = "*kʷeys-"
    pie_causative = "*kʷoyséye-"
    proto_germanic = "*hʷayzjaną"
    old_english = "hwǣreþ"
    middle_english = "whereth"

    print("This script traces a hypothetical word from Proto-Indo-European to Middle English.")
    print("-" * 70)

    # Step 1: Proto-Indo-European to Proto-Germanic
    print("1. From Proto-Indo-European (PIE) to Proto-Germanic (PGmc):")
    print(f"   - Start with PIE root: {pie_root} ('to see, to heed').")
    print(f"   - Form the o-grade causative: {pie_causative} ('to make see' -> 'to show').")
    print("   - Apply sound changes:")
    print("     - Grimm's Law: *kʷ > *hʷ")
    print("     - Verner's Law: *s > *z (due to accent on the suffix)")
    print("     - Vowel/Suffix Shift: *oy > *ay; *-éye- > *-janą (Class 1 Weak Verb)")
    print(f"   - Resulting PGmc infinitive: {proto_germanic}\n")

    # Step 2: Proto-Germanic to Old English
    print("2. From Proto-Germanic (PGmc) to Old English (OE):")
    print(f"   - Start with PGmc infinitive: {proto_germanic}.")
    print("   - Apply sound changes:")
    print("     - I-umlaut: The /j/ in *-janą fronts the root diphthong *ai > ǣ.")
    print("     - Rhotacism: Medial *z > r.")
    print("     - The infinitive ending becomes -an, giving 'hwǣran'.")
    print("   - For the 3rd person singular present, add the suffix '-eþ'.")
    print(f"   - Resulting Old English form: {old_english}\n")

    # Step 3: Old English to Middle English
    print("3. From Old English (OE) to Middle English (ME):")
    print(f"   - Start with OE form: {old_english}.")
    print("   - Apply sound changes:")
    print("     - Vowel Shift: Long ǣ > long ē (/e:/).")
    print("     - Spelling Change: hw > wh.")
    print("     - Inflectional Ending: The ending -eþ becomes -eth.")
    print(f"   - Resulting Middle English form: {middle_english}\n")
    
    print("-" * 70)
    print("The final derived third person singular present form in Middle English, meaning 'he shows', would be:")
    print(f"{middle_english}")

trace_word_evolution()
<<<whereth>>>