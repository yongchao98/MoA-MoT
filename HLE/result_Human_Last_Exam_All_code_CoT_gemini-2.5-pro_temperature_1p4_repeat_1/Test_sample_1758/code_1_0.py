def transcribe_to_puerto_rican_spanish():
    """
    This script explains and constructs the phonetic transcription for
    the English borrowing "shock absorber" in Puerto Rican Spanish.
    """
    print("Step 1: The original English word is 'shock absorber', pronounced roughly as /ʃɑk əbˈsɔɹbɚ/.")
    print("Step 2: When borrowed into Spanish, it's commonly shortened and blended into a single word like 'chosorber'.")
    print("\nApplying the phonological rules of Puerto Rican Spanish to 'chosorber':\n")

    # Breaking down the word into its phonetically relevant parts for transformation
    part1 = "'cho-'"
    part2 = "'-s-'"
    part3 = "'-or-'"
    part4 = "'-ber'"

    # Applying rules to each part
    result1_char = "tʃo"
    result1_explanation = "(The 'sh' sound /ʃ/ becomes 'ch' /tʃ/, and the vowel 'o' is used)."

    result2_char = "h"
    result2_explanation = "(Syllable-final 's' is aspirated, sounding like an 'h')."

    result3_char = "ol"
    result3_explanation = "(The 'r' in 'or' is at the end of a syllable, so it is lateralized to an 'l')."

    result4_char = "ˈβel"
    result4_explanation = "(The 'b' is between vowels, so it softens to [β]. The final 'r' is also lateralized to 'l'. The stress falls on this last syllable)."

    print(f"1. The 'shock' part becomes '{result1_char}': {result1_explanation}")
    print(f"2. The 's' from 'absorber' becomes '{result2_char}': {result2_explanation}")
    print(f"3. The 'or' from 'absorber' becomes '{result3_char}': {result3_explanation}")
    print(f"4. The final 'ber' becomes '{result4_char}': {result4_explanation}")
    
    # Constructing the final phonetic transcription
    final_transcription_symbols = ["tʃ", "o", "h", ".", "o", "l", "ˈ", "β", "e", "l"]
    final_transcription_string = "[tʃoh.olˈβel]"
    
    print("\nPutting it all together, the final phonetic transcription is assembled from these parts.")
    print("Final Pronunciation: " + final_transcription_string)
    
    # Printing each symbol of the final transcription as requested
    print("\nThe symbols in the final transcription are:")
    for symbol in final_transcription_symbols:
        print(f"'{symbol}'")

transcribe_to_puerto_rican_spanish()

print("\n<<<[tʃoh.olˈβel]>>>")