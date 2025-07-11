def find_similar_words():
    """
    This function presents the answer by printing the comparison between
    two languages from different Asian cultures: Malay and Gujarati.
    """

    language1 = "Malay"
    culture1 = "Austronesian (spoken in Malaysia, Indonesia)"
    language2 = "Gujarati"
    culture2 = "Indo-Aryan (spoken in Gujarat, India)"

    print(f"Finding two languages in Asia with no direct contact but similar words for 'dad', 'mom', and 'broom'.")
    print("="*80)
    print(f"The chosen languages are {language1} and {language2}.")
    print("\n[CULTURAL & LINGUISTIC BACKGROUND]")
    print(f"{language1} is an Austronesian language, while {language2} is an Indo-Aryan language.")
    print("They belong to different language families and their speakers are geographically distant with no direct historical contact, making the similarities remarkable.")

    # Define the word data in a dictionary
    word_comparisons = {
        'Dad': {
            'Language 1': f"{language1}",
            'Word 1': "bapa",
            'Language 2': f"{language2}",
            'Word 2': "bāpuji",
            'Note': "Both share the distinct 'bap-' root for father."
        },
        'Mom': {
            'Language 1': f"{language1}",
            'Word 1': "emak",
            'Language 2': f"{language2}",
            'Word 2': "mā",
            'Note': "Both are derived from the universal 'ma' nursery sound for mother."
        },
        'Broom': {
            'Language 1': f"{language1}",
            'Word 1': "penyapu",
            'Language 2': f"{language2}",
            'Word 2': "saavarni",
            'Note': "The root 'sapu' (to sweep) in Malay is phonetically similar to the root 'saav-' in Gujarati."
        }
    }

    print("\n[WORD COMPARISON]")
    # Print the comparisons in a structured way
    for term, data in word_comparisons.items():
        print(f"\n--- Word: {term} ---")
        print(f"{data['Language 1']:<10}: {data['Word 1']:<15}")
        print(f"{data['Language 2']:<10}: {data['Word 2']:<15}")
        print(f"Note: {data['Note']}")

if __name__ == "__main__":
    find_similar_words()
