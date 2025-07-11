def find_similar_words():
    """
    Identifies and prints a comparison of words from two Asian languages
    with no direct historical contact.
    """
    language1 = "Korean"
    language2 = "Persian (Farsi)"
    
    word_comparisons = {
        "mom": {
            "lang1_word": "eomma (엄마)",
            "lang2_word": "mâmân (مامان)"
        },
        "dad": {
            "lang1_word": "appa (아빠)",
            "lang2_word": "bābā (بابا)"
        },
        "broom": {
            "lang1_word": "bitjaru (빗자루)",
            "lang2_word": "jâru (جارو)"
        }
    }

    print(f"Two languages from two different Asian cultures with no direct contact are {language1} and {language2}.")
    print("One is from East Asia and the other from West Asia, and they belong to different language families.")
    print("\nBelow is a comparison of the words for 'mom', 'dad', and 'broom':\n")
    
    # Print each comparison in a format resembling an equation
    for english_word, details in word_comparisons.items():
        korean_word = details["lang1_word"]
        persian_word = details["lang2_word"]
        
        print(f"English word: '{english_word}'")
        print("--- Comparison ---")
        print(f"{language1}: {korean_word}")
        print(f"{language2}: {persian_word}")
        print("==================\n")

if __name__ == "__main__":
    find_similar_words()