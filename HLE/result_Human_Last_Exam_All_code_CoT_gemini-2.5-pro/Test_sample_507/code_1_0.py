def find_similar_words():
    """
    This function identifies and prints information about two languages from Asia
    with similar words for "mom", "dad", and "broom" despite having no
    direct cultural contact.
    """
    
    print("The two languages are Korean and Farsi (Persian).")
    print("These languages originate from very different cultures in East Asia and West Asia, respectively, with no direct historical contact between them.\n")

    # Define the words for each language
    lang1_name = "Korean"
    lang1_mom = "엄마 (eomma)"
    lang1_dad = "아빠 (appa)"
    lang1_broom = "빗자루 (bitjaru)"

    lang2_name = "Farsi (Persian)"
    lang2_mom = "مامان (mâmân)"
    lang2_dad = "بابا (bâbâ)"
    lang2_broom = "جارو (jâru)"

    # Print the comparison in a clear format
    print("--- Word Comparison ---\n")

    # Mom
    print("Category: Mom")
    print(f"{lang1_name}: {lang1_mom}")
    print(f"{lang2_name}:  {lang2_mom}")
    print("Analysis: Both words are variations of the universal 'mama' nursery word.\n")

    # Dad
    print("Category: Dad")
    print(f"{lang1_name}: {lang1_dad}")
    print(f"{lang2_name}:  {lang2_dad}")
    print("Analysis: Both words are variations of the universal 'papa'/'baba' nursery word.\n")

    # Broom
    print("Category: Broom")
    print(f"{lang1_name}: {lang1_broom}")
    print(f"{lang2_name}:  {lang2_broom}")
    print("Analysis: This is a remarkable coincidence. The Korean word is a compound of 'bit' (comb) + 'jaru' (handle). The Farsi word 'jâru' is phonetically almost identical to the 'jaru' component in the Korean word.\n")

    # Final summary as requested, printing each word again.
    print("--- Final Equation of Similarity ---")
    print(f"Words in {lang1_name}: '{lang1_mom}', '{lang1_dad}', '{lang1_broom}'")
    print("are highly similar to")
    print(f"Words in {lang2_name}: '{lang2_mom}', '{lang2_dad}', '{lang2_broom}'")

if __name__ == '__main__':
    find_similar_words()
