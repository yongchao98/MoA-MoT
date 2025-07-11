def identify_languages():
    """
    This function identifies two languages based on a set of linguistic clues
    and prints the solution.
    """

    # --- Language A ---
    lang_a = "Scottish Gaelic"
    reasoning_a = [
        "1. No 'k' or 'w': The traditional Scottish Gaelic alphabet has 18 letters and does not include 'k' or 'w'.",
        "2. Has 'à': The grave accent is used on all five vowels in Gaelic, including 'à'.",
        "3. Recognized minority language: It is an officially recognized minority language in the UK (Scotland)."
    ]

    # --- Language B ---
    lang_b = "Northern Sami"
    reasoning_b = [
        "1. Widely used 'ggj' and 'skt': These consonant clusters are common in Northern Sami due to its system of consonant gradation.",
        "   - 'ggj' appears in words like 'loggje' (to lay).",
        "   - 'skt' appears in words like 'básket' (to wash).",
        "2. Recognized minority language: It is an officially recognized minority language in Norway, Sweden, and Finland."
    ]

    # --- Print the results ---
    print("Based on the clues provided:")
    print("-" * 30)

    print(f"Language a is: {lang_a}")
    print("Reasoning:")
    for reason in reasoning_a:
        print(f"  - {reason}")

    print("\n" + "-" * 30)

    print(f"Language b is: {lang_b}")
    print("Reasoning:")
    for reason in reasoning_b:
        print(f"  - {reason}")
    print("-" * 30)

if __name__ == "__main__":
    identify_languages()