def find_linguistic_coincidence():
    """
    This script presents a famous linguistic coincidence between Hungarian and
    Mongolian, two languages from cultures with Asian roots but no direct
    linguistic ancestry, which share surprisingly similar words for
    "mom", "dad", and "broom".
    """

    print("Comparing similar core vocabulary in Hungarian and Mongolian:\n")

    # Data structure: (English Word, Hungarian Word, Mongolian Word)
    word_pairs = [
        ("Mom", "Anya", "Eej (pronounced similarly to 'Anj')"),
        ("Dad", "Apa", "Aav"),
        ("Broom", "Seprű", "Shüür")
    ]

    for english, hungarian, mongolian in word_pairs:
        print(f"The word for \"{english}\":")
        print(f"  - In Hungarian is: {hungarian}")
        print(f"  - In Mongolian is: {mongolian}")
        print("-" * 20)

find_linguistic_coincidence()