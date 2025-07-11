def find_similar_words():
    """
    This script identifies two languages from different Asian cultures with no direct contact
    that share similar words for "mom," "dad," and "broom."
    """

    lang1 = "Malay"
    lang2 = "Turkish"
    culture1 = "Maritime Southeast Asia"
    culture2 = "Anatolia / West Asia"

    print(f"The two languages are {lang1} and {lang2}.")
    print(f"These languages are from two distant Asian cultures ({culture1} and {culture2})")
    print("with unrelated linguistic origins (Austronesian and Turkic) and no significant direct historical contact.")
    print("-" * 20)
    print("\nHere is the word comparison:\n")

    # Mom
    mom1 = "Mama"
    mom2 = "Mama"
    print(f"1. Word for 'Mom':")
    print(f"   In both {lang1} and {lang2}, the colloquial nursery word for mother is nearly identical.")
    print(f"   - {lang1}: '{mom1}' (colloquial/nursery term)")
    print(f"   - {lang2}: '{mom2}' (colloquial/nursery term)")
    print("")

    # Dad
    dad1 = "Bapa"
    dad2 = "Baba"
    print(f"2. Word for 'Dad':")
    print(f"   The word for father is also highly similar.")
    print(f"   - {lang1}: '{dad1}'")
    print(f"   - {lang2}: '{dad2}'")
    print("")

    # Broom
    broom1_noun = "penyapu"
    broom1_verb = "sapu"
    broom2_noun = "süpürge"
    broom2_verb = "süpürmek"
    print(f"3. Word for 'Broom':")
    print(f"   The nouns for 'broom' ({broom1_noun} vs. {broom2_noun}) are different.")
    print(f"   However, the similarity is high in the verb root for 'to sweep', from which the word for 'broom' is derived.")
    print(f"   - {lang1} root for 'sweep': '{broom1_verb}'")
    print(f"   - {lang2} root for 'sweep': '{broom2_verb}' (root: süpür-)")
    print("")

    print("-" * 20)
    print("Final comparison:")
    print(f"{lang1} 'Mom' ({mom1}) == {lang2} 'Mom' ({mom2})")
    print(f"{lang1} 'Dad' ({dad1}) ~= {lang2} 'Dad' ({dad2})")
    print(f"{lang1} 'Broom' (root: {broom1_verb}) ~= {lang2} 'Broom' (root: süpür-)")


if __name__ == "__main__":
    find_similar_words()

<<<Malay and Turkish>>>