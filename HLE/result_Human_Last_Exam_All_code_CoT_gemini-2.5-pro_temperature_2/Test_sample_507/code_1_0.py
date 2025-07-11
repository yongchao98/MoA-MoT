def find_similar_words():
    """
    Identifies and prints similar words from two distant Asian languages: Korean and Tamil.
    These languages belong to unrelated language families (Koreanic and Dravidian)
    and are geographically separated, making the similarities in core family vocabulary remarkable.
    """
    # Data structure: {English Word: (Language1, Romanization1, Language2, Romanization2)}
    word_comparisons = {
        "Mom": ("엄마", "eomma", "அம்மா", "ammā"),
        "Dad": ("아빠", "appa", "அப்பா", "appā"),
        "Broom": ("빗자루", "bitjaru", "விளக்குமாறு", "vilakkumaru")
    }

    print("Comparing words in Korean and Tamil:\n")
    print(f"{'English':<10} | {'Korean':<20} | {'Tamil':<20}")
    print("-" * 55)

    for word, (korean_word, korean_roman, tamil_word, tamil_roman) in word_comparisons.items():
        korean_display = f"{korean_word} ({korean_roman})"
        tamil_display = f"{tamil_word} ({tamil_roman})"
        print(f"{word:<10} | {korean_display:<20} | {tamil_display:<20}")

    print("\nNote: The words for 'Mom' and 'Dad' are nearly identical.")
    print("The word for 'Broom' is not similar and is included to fulfill the prompt's requirements,")
    print("highlighting that the coincidence is specific to certain words.")


if __name__ == "__main__":
    find_similar_words()