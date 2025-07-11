def find_similar_words():
    """
    This function identifies and prints a comparison of words for "mom", "dad",
    and "broom" from two geographically and culturally distant Asian languages.
    """

    languages = {
        "Korean": {
            "culture": "East Asian",
            "language_family": "Koreanic",
            "mom": "eomma (엄마)",
            "dad": "appa (아빠)",
            "broom": "bitjaru (빗자루)"
        },
        "Tamil": {
            "culture": "South Asian (Dravidian)",
            "language_family": "Dravidian",
            "mom": "ammā (அம்மா)",
            "dad": "appā (அப்பா)",
            "broom": "thudaippam (துடைப்பம்)"
        }
    }

    lang1_name = "Korean"
    lang2_name = "Tamil"
    lang1 = languages[lang1_name]
    lang2 = languages[lang2_name]

    print(f"Finding languages with similar words for 'mom', 'dad', and 'broom' with no direct contact:\n")
    print(f"Suggested Languages: {lang1_name} and {lang2_name}")
    print(f"Rationale: These languages are from different families ({lang1['language_family']} vs. {lang2['language_family']}) and cultures ({lang1['culture']} vs. {lang2['culture']}) with no significant direct contact.")
    print("-" * 50)

    print(f"Comparing the words:\n")

    # Mom comparison
    print("Word: Mom")
    print(f"  {lang1_name}: {lang1['mom']}")
    print(f"  {lang2_name}: {lang2['mom']}")
    print(f"  Similarity: Very High. Both use a similar vowel sound and a double 'm'.\n")

    # Dad comparison
    print("Word: Dad")
    print(f"  {lang1_name}: {lang1['dad']}")
    print(f"  {lang2_name}: {lang2['dad']}")
    print(f"  Similarity: Very High. Both use a similar vowel sound and a double 'p'.\n")

    # Broom comparison
    print("Word: Broom")
    print(f"  {lang1_name}: {lang1['broom']}")
    print(f"  {lang2_name}: {lang2['broom']}")
    print(f"  Similarity: Low. These words are not similar, indicating the coincidence may be limited to parental terms.")
    print("-" * 50)


if __name__ == '__main__':
    find_similar_words()