def find_similar_words():
    """
    This script identifies two languages from different Asian cultures
    with similar words for 'mom', 'dad', and 'broom'.
    The chosen languages are Malay and Telugu, which are from unrelated
    language families and have no significant history of direct contact.
    """

    # Data for the two languages
    # Language 1: Malay (Austronesian family)
    malay_words = {
        "language": "Malay",
        "culture": "Malay Archipelago (Southeast Asia)",
        "mom": "emak",
        "dad": "bapa",
        "broom": "penyapu"
    }

    # Language 2: Telugu (Dravidian family)
    # 'appa' is a common word for father, alongside 'nanna'.
    telugu_words = {
        "language": "Telugu",
        "culture": "South-Central India (South Asia)",
        "mom": "amma",
        "dad": "appa",
        "broom": "cheepuru"
    }

    print(f"Comparing two languages from Asia with no direct contact: {malay_words['language']} and {telugu_words['language']}.\n")
    print(f"{malay_words['language']} is an Austronesian language from {malay_words['culture']}.")
    print(f"{telugu_words['language']} is a Dravidian language from {telugu_words['culture']}.\n")
    print("The following similarities are coincidental, not from a shared origin.\n")
    print("--- Comparison ---")

    # We will represent the comparisons in a simple format.
    # The numbers '1' and '2' simply pair the words for comparison.
    print(f"1. For Mom: '{malay_words['mom']}' (Malay) is similar to '{telugu_words['mom']}' (Telugu).")
    print(f"2. For Dad: '{malay_words['dad']}' (Malay) is similar to '{telugu_words['dad']}' (Telugu).")
    print(f"3. For Broom: '{malay_words['broom']}' (Malay) is similar to '{telugu_words['broom']}' (Telugu).")

find_similar_words()