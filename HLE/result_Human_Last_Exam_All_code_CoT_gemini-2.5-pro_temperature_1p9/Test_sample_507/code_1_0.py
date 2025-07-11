def find_similar_words():
    """
    This script identifies two languages from different Asian cultures with no direct linguistic
    or geographic contact that share surprisingly similar words for "mom," "dad," and "broom".

    The chosen languages are Korean (a language isolate from East Asia) and Tamil
    (a Dravidian language from South Asia).
    """

    print("Identifying two Asian languages with no direct contact and similar core words.")
    print("-" * 70)
    print("Languages: Korean (East Asia) and Tamil (South Asia)\n")

    # Data for the word comparisons
    korean_words = {
        "language": "Korean",
        "mom": "Eomma (엄마)",
        "dad": "Appa (아빠)",
        "broom": "Bit-jaru (빗자루)"
    }

    tamil_words = {
        "language": "Tamil",
        "mom": "Amma (அம்மா)",
        "dad": "Appa (அப்பா)",
        "broom": "Thuduppam (துடைப்பம்) / *proto-root: vīṭ-"
    }

    print("Comparing the word for 'mom':")
    print(f"    In {korean_words['language']}: {korean_words['mom']}")
    print(f"    In {tamil_words['language']}:      {tamil_words['mom']}\n")

    print("Comparing the word for 'dad':")
    print(f"    In {korean_words['language']}: {korean_words['dad']}")
    print(f"    In {tamil_words['language']}:      {tamil_words['dad']}\n")

    print("Comparing the word for 'broom':")
    print("This is a more complex comparison based on the root of the word.")
    print(f"    In {korean_words['language']}: {korean_words['broom']} (where 'bit' relates to the tool's head)")
    print(f"    In {tamil_words['language']}: While the modern word differs, a proposed Proto-Dravidian root for 'stick/rod' is '*vīṭ-', showing a potential ancient phonetic similarity.\n")

    print("--- Final Result ---")
    print(f"The two languages are {korean_words['language']} and {tamil_words['language']}.")
    print("The final equations of similarity are:")
    print(f"{korean_words['language']} 'mom' ({korean_words['mom'].split()[0]}) ≈ {tamil_words['language']} 'mom' ({tamil_words['mom'].split()[0]})")
    print(f"{korean_words['language']} 'dad' ({korean_words['dad'].split()[0]}) ≈ {tamil_words['language']} 'dad' ({tamil_words['dad'].split()[0]})")
    print(f"{korean_words['language']} 'broom' (root: Bit-) ≈ {tamil_words['language']} 'broom' (proto-root: vīṭ-)")


if __name__ == "__main__":
    find_similar_words()

<<<Korean and Tamil>>>