def find_similar_language_words():
    """
    Identifies and prints two Asian languages with no direct contact
    that have surprisingly similar words for "broom", "mom", and "dad".
    """

    language1 = {
        "name": "Korean",
        "culture_location": "East Asia",
        "words": {
            "mom": "eomma (엄마)",
            "dad": "appa (아빠)",
            "broom": "bitjaru (빗자루)"
        }
    }

    language2 = {
        "name": "Persian (Farsi)",
        "culture_location": "West Asia",
        "words": {
            "mom": "māmān (مامان)",
            "dad": "bābā (بابا)",
            "broom": "jāru (جارو)"
        }
    }

    print(f"Finding two languages from different Asian cultures with no direct contact and similar core words.\n")
    print("="*70)
    print(f"Selected Languages: {language1['name']} and {language2['name']}\n")
    print(f"Analysis of Word Similarities:")
    print("="*70)

    # Mom
    print(f"1. Word for 'Mom':")
    print(f"   - Korean: {language1['words']['mom']}")
    print(f"   - Persian: {language2['words']['mom']}")
    print("   - Similarity: Both are nursery words based on the simple, universal 'ma' sound.\n")

    # Dad
    print(f"2. Word for 'Dad':")
    print(f"   - Korean: {language1['words']['dad']}")
    print(f"   - Persian: {language2['words']['dad']}")
    print("   - Similarity: Both use a repeated bilabial consonant ('p'/'b') and the 'a' vowel, a common pattern in nursery words for 'father'.\n")

    # Broom
    print(f"3. Word for 'Broom':")
    print(f"   - Korean: {language1['words']['broom']}")
    print(f"   - Persian: {language2['words']['broom']}")
    print("   - Similarity: This is the most striking match. The Korean word 'bitjaru' contains a segment, 'jaru', that is nearly identical to the complete Persian word 'jāru'.")
    print("="*70)

if __name__ == '__main__':
    find_similar_language_words()