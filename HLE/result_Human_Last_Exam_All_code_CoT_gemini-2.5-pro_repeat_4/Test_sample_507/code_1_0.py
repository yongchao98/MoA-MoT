def find_similar_words():
    """
    This function presents a solution to the linguistic riddle by identifying
    two languages from different Asian cultures with similar words for "mom",
    "dad", and "broom".
    """

    language1 = "Malay"
    language2 = "Turkish"

    dad_malay = "Bapa"
    dad_turkish = "Baba"

    mom_malay = "Emak"
    mom_turkish = "Ana"

    broom_malay = "Penyapu"
    broom_turkish = "Süpürge"

    verb_malay = "sapu"
    verb_turkish = "süpür"

    print(f"Two languages from two cultures in Asia with no direct contact that share highly similar words are {language1} and {language2}.")
    print("-" * 80)
    print("1. The Languages and Cultures:")
    print(f"   - {language1}: An Austronesian language from Southeast Asia.")
    print(f"   - {language2}: A Turkic language from West Asia.")
    print("   These cultures are geographically distant and belong to different language families, with minimal historical contact.")
    print("-" * 80)
    print("2. The Word Comparisons:\n")

    print("   Word: DAD")
    print(f"   - In {language1}: {dad_malay}")
    print(f"   - In {language2}: {dad_turkish}")
    print("   -> These two words are nearly identical.\n")

    print("   Word: MOM")
    print(f"   - In {language1} (common, informal): {mom_malay}")
    print(f"   - In {language2} (word for 'mother'): {mom_turkish}")
    print("   -> These words share similar phonetic structures common for 'mother' words.\n")

    print("   Word: BROOM")
    print(f"   - In {language1}: {broom_malay}")
    print(f"   - In {language2}: {broom_turkish}")
    print("   -> The similarity is found in the root verbs for 'to sweep', from which 'broom' is derived:")
    print(f"      - {language1} root verb: '{verb_malay}'")
    print(f"      - {language2} root verb: '{verb_turkish}'")
    print("   -> The verbs are phonetically very similar, making this a remarkable coincidence.")

find_similar_words()