def present_linguistic_comparison():
    """
    This script presents a comparison between Korean and Tamil,
    two unrelated Asian languages, to show surprising similarities
    in core vocabulary, as requested.
    """

    # --- Language Data ---
    lang1 = {
        "name": "Korean",
        "culture": "East Asian",
        "family": "Koreanic",
        "mom": "eomma (엄마)",
        "dad": "appa (아빠)",
        "broom": "bitjaru (빗자루)"
    }

    lang2 = {
        "name": "Tamil",
        "culture": "South Asian (Dravidian)",
        "family": "Dravidian",
        "mom": "ammā (அம்மா)",
        "dad": "appā (அப்பா)",
        "broom": "thudaippam (துடைப்பம்)"
    }

    # --- Presentation ---
    print("Language Comparison: Korean vs. Tamil\n")
    print(f"This comparison shows two languages from different Asian cultures with no direct contact.")
    print("-" * 40)
    print(f"Language 1: {lang1['name']} ({lang1['culture']}, {lang1['family']} family)")
    print(f"Language 2: {lang2['name']} ({lang2['culture']}, {lang2['family']} family)")
    print("-" * 40 + "\n")

    print("Comparing the word for 'mom':")
    print(f"{lang1['name']}: {lang1['mom']}")
    print(f"{lang2['name']}: {lang2['mom']}")
    print("Result: The words are highly similar.\n")

    print("Comparing the word for 'dad':")
    print(f"{lang1['name']}: {lang1['dad']}")
    print(f"{lang2['name']}: {lang2['dad']}")
    print("Result: The words are nearly identical.\n")

    print("Comparing the word for 'broom':")
    print(f"{lang1['name']}: {lang1['broom']}")
    print(f"{lang2['name']}: {lang2['broom']}")
    print("Result: The words are not similar, as is expected for unrelated languages.")

if __name__ == '__main__':
    present_linguistic_comparison()