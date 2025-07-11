import unicodedata

def generate_transcription():
    """
    This function explains and generates the phonetic transcription of "shock absorber"
    in Puerto Rican Spanish.
    """

    # Define IPA characters for clarity
    SH = unicodedata.lookup("LATIN SMALL LETTER ESH")
    CH = "t" + unicodedata.lookup("COMBINING TILDE") # A way to represent t͡ʃ, but simpler is 'tʃ'
    CH_AFFRICATE = "tʃ"
    STRESS_MARK = "ˈ"
    SCHWA = "ə"
    OPEN_O = "ɔ"
    ALPHA = "ɑ"

    print("This script derives the Puerto Rican Spanish pronunciation of the English phrase 'shock absorber'.")
    print("------------------------------------------------------------------------------------------------\n")

    # Step 1: Adaptation of "shock"
    print("Part 1: Adapting 'shock' (English IPA: /ˈʃɑk/)")
    print(f"1. The English 'sh' sound [{SH}] is adapted to the Spanish 'ch' sound [{CH_AFFRICATE}].")
    print(f"2. The English vowel [{ALPHA}] (as in 'shock') is pronounced as the Spanish vowel [o].")
    shock_result = f"[{CH_AFFRICATE}ok]"
    print(f"--> The adapted pronunciation of 'shock' is: {shock_result}\n")

    # Step 2: Adaptation of "absorber"
    print("Part 2: Adapting 'absorber' (English IPA: /əbˈsɔrbər/)")
    print("1. Apheresis: The initial unstressed syllable 'ab-' is dropped, a common feature in borrowings. The word is simplified to 'sorber'.")
    print(f"2. Lateralization: The syllable-final /r/ sounds are pronounced as [l]. This is a key feature of Puerto Rican Spanish.")
    print(f"   - 'sorb' -> [solb]")
    print(f"   - '-er' -> [el]")
    print(f"3. Vowel Mapping: English vowel [{OPEN_O}] becomes [o] and [{SCHWA}] becomes [e].")
    absorber_result = "[solbel]"
    print(f"--> The adapted pronunciation of 'sorber' is: {absorber_result}\n")

    # Step 3: Final Combination
    print("Part 3: Final Combined Transcription")
    print("The two parts are combined, retaining stress on the original stressed syllables.")
    final_transcription = f"[{STRESS_MARK}{CH_AFFRICATE}ok {STRESS_MARK}solbel]"
    print(f"\nFinal phonetic transcription: {final_transcription}\n")

    print("Breaking down the final transcription into its phonemic components:")
    # The prompt requests that each "number" (interpreted as "phoneme") be output.
    final_phonemes = [
        STRESS_MARK + CH_AFFRICATE, 'o', 'k',
        ' ',
        STRESS_MARK + 's', 'o', 'l', 'b', 'e', 'l'
    ]
    # Adding brackets for the full notation
    print("[", end=" ")
    for phoneme in final_phonemes:
        print(phoneme, end=" ")
    print("]")


# Execute the function to get the answer.
generate_transcription()
