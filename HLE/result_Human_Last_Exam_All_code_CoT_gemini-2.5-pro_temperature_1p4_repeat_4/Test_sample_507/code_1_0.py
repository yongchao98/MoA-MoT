def find_similar_words():
    """
    This script identifies two languages from different Asian cultures with no direct contact
    where the words for "mom," "dad," and "broom" show surprising similarities.
    """

    # Language 1: Korean
    lang1 = "Korean"
    culture1 = "East Asia (Korean Peninsula)"
    lang1_mom = "Eomma (엄마)"
    lang1_dad = "Appa (아빠)"
    lang1_broom = "Bitjaru (빗자루)"

    # Language 2: Telugu (a Dravidian language)
    lang2 = "Telugu"
    culture2 = "South-Central Asia (Indian Subcontinent)"
    lang2_mom = "Amma (అమ్మ)"
    lang2_dad = "Appa (అప్ప)"
    lang2_broom = "Cheepuru (చీపురు)"

    print(f"Finding two languages from different Asian cultures with similar words for 'mom,' 'dad,' and 'broom'.\n")
    print(f"The chosen languages are {lang1} and {lang2}, which are from different language families and have no direct historical contact.\n")

    print("-" * 40)
    print(f"{'Word':<10} | {'Language 1: Korean':<25} | {'Language 2: Telugu':<25}")
    print("-" * 40)
    print(f"{'Mom':<10} | {lang1_mom:<25} | {lang2_mom:<25}")
    print(f"{'Dad':<10} | {lang1_dad:<25} | {lang2_dad:<25}")
    print(f"{'Broom':<10} | {lang1_broom:<25} | {lang2_broom:<25}")
    print("-" * 40)
    print("\nObservation:")
    print(f"- 'Mom': {lang1_mom} and {lang2_mom} are nearly identical.")
    print(f"- 'Dad': {lang1_dad} and {lang2_dad} are also nearly identical.")
    print(f"- 'Broom': {lang1_broom} and {lang2_broom} are not identical but share some phonetic characteristics, such as the 'ru' sound at the end and a 'p'/'b' consonant.")

# Execute the function to print the result
find_similar_words()
