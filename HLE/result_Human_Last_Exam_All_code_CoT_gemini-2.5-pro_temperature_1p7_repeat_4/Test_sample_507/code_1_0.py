def find_similar_words_in_asian_languages():
    """
    This function identifies and prints two languages from different Asian cultures
    with no direct contact that share similar words for "broom", "mom", and "dad".
    """
    language1 = "Hindi"
    culture1 = "Indian (Indo-Aryan)"
    language2 = "Korean"
    culture2 = "Korean (Koreanic)"

    print(f"The two languages are {language1} and {language2}.")
    print("-" * 30)
    print(f"Cultures and Contact: The {culture1} culture of South Asia and the {culture2} culture of East Asia are geographically distant and belong to completely different language families. They have no history of significant direct contact that would explain borrowing of such basic vocabulary.")
    print("-" * 30)
    print("Here is the word comparison:\n")

    # Word data
    words = {
        "Broom": {"Hindi": "jhāṛū (झाड़ू)", "Korean": "bitjaru (빗자루)"},
        "Mom": {"Hindi": "mā (माँ)", "Korean": "eomma (엄마)"},
        "Dad": {"Hindi": "bāp (बाप)", "Korean": "appa (아빠)"}
    }

    print(f"{'Word':<10} | {'Hindi':<15} | {'Korean':<15} | {'Similarity Analysis'}")
    print(f"{'-'*10} | {'-'*15} | {'-'*15} | {'-'*30}")

    # Broom
    term = "Broom"
    hindi_word = words[term]["Hindi"]
    korean_word = words[term]["Korean"]
    analysis_broom = "The Hindi 'jhāṛū' and the root 'jaru' in Korean 'bitjaru' are phonetically very similar."
    print(f"{term:<10} | {hindi_word:<15} | {korean_word:<15} | {analysis_broom}")

    # Mom
    term = "Mom"
    hindi_word = words[term]["Hindi"]
    korean_word = words[term]["Korean"]
    analysis_mom = "Both 'mā' and 'eomma' are archetypal nursery words based on the 'ma' sound."
    print(f"{term:<10} | {hindi_word:<15} | {korean_word:<15} | {analysis_mom}")

    # Dad
    term = "Dad"
    hindi_word = words[term]["Hindi"]
    korean_word = words[term]["Korean"]
    analysis_dad = "Both 'bāp' and 'appa' are nursery words using similar bilabial sounds (b/p)."
    print(f"{term:<10} | {hindi_word:<15} | {korean_word:<15} | {analysis_dad}")

if __name__ == "__main__":
    find_similar_words_in_asian_languages()