def find_similar_words_in_asian_languages():
    """
    This script identifies and explains the similarity in the words for "mom",
    "dad", and "broom" between Hindi and Vietnamese, two languages from
    different Asian cultures with no direct historical contact.
    """

    print("Finding two languages from different Asian cultures with no direct contact where the words for 'mom', 'dad', and 'broom' are highly similar.")
    print("="*80)

    lang1 = {
        "name": "Hindi",
        "culture": "Indian (South Asia)",
        "family": "Indo-Aryan",
        "words": {
            "mom": "Mā (माँ)",
            "dad": "Bāp (बाप)",
            "broom": "Jhāṛū (झाड़ू)"
        }
    }

    lang2 = {
        "name": "Vietnamese",
        "culture": "Vietnamese (Southeast Asia)",
        "family": "Austroasiatic",
        "words": {
            "mom": "Mẹ",
            "dad": "Ba",
            "broom": "Chổi"
        }
    }

    print(f"Language 1: {lang1['name']} ({lang1['culture']})")
    print(f"Language 2: {lang2['name']} ({lang2['culture']})")
    print("-" * 80)

    print("Word Comparison:\n")

    # The prompt requested to output each part of the "equation",
    # so we will output each word pair clearly.
    print("1. The word for 'mom':")
    print(f"   - In {lang1['name']}: '{lang1['words']['mom']}'")
    print(f"   - In {lang2['name']}: '{lang2['words']['mom']}'")
    print("   -> Similarity: Both are short, nasal words starting with 'M'.\n")

    print("2. The word for 'dad':")
    print(f"   - In {lang1['name']}: '{lang1['words']['dad']}'")
    print(f"   - In {lang2['name']}: '{lang2['words']['dad']}'")
    print("   -> Similarity: Both are simple words starting with a 'B' sound.\n")

    print("3. The word for 'broom':")
    print(f"   - In {lang1['name']}: '{lang1['words']['broom']}'")
    print(f"   - In {lang2['name']}: '{lang2['words']['broom']}'")
    print("   -> Similarity: Both start with a similar-sounding consonant ([d͡ʒʱ] vs [c]).\n")

    print("-" * 80)
    print("Conclusion:")
    print(f"{lang1['name']} and {lang2['name']} fit the criteria. They belong to different language families and evolved in different regions of Asia with no direct contact that would explain borrowing these words.")
    print("The similarities are due to:")
    print("  - 'Mom' and 'Dad': A global linguistic phenomenon known as 'nursery words'.")
    print("  - 'Broom': Pure coincidence.")

if __name__ == '__main__':
    find_similar_words_in_asian_languages()