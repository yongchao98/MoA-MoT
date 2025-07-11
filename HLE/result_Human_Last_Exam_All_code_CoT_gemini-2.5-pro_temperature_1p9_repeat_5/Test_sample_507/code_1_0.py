def find_similar_words():
    """
    This function identifies and prints information about two languages from two
    cultures in Asia with no direct contact, where the words for "broom", "mom",
    and "dad" are nearly identical or highly similar.
    """

    print("The two languages selected are Persian (Farsi) and Malagasy.")
    print("-" * 50)
    print("1. Language & Culture: Persian (Farsi)")
    print("   - Region: West Asia (Iran)")
    print("   - Language Family: Indo-European")
    print("\n2. Language & Culture: Malagasy")
    print("   - Region: Madagascar. The culture's origins are from Austronesian peoples of Southeast Asia (Borneo).")
    print("   - Language Family: Austronesian")
    print("\n   - Condition Check: These two cultures are geographically distant and have no significant history of direct contact.")
    print("-" * 50)
    print("Comparing the words for 'broom', 'mom', and 'dad':\n")

    # Word data in a list of dictionaries for easy iteration
    word_data = [
        {"concept": "Broom", "persian": "jâru", "malagasy": "tsarò"},
        {"concept": "Mom", "persian": "mâmân", "malagasy": "mama"},
        {"concept": "Dad", "persian": "bâbâ", "malagasy": "dada"}
    ]

    # The prompt "output each number in the final equation!" is interpreted as
    # "output each word in the final comparison". The following loop does this.
    for item in word_data:
        print(f"For the concept '{item['concept']}', the words are:")
        # This line prints the comparison like a pseudo-equation
        print(f"  Persian '{item['persian']}' ≈ Malagasy '{item['malagasy']}'\n")

    # Construct the final answer string for the required format
    final_answer = ("The two languages are Persian (Farsi) from West Asia and Malagasy (with cultural roots in Southeast Asia). "
                    "The words are: Broom (Persian: jâru, Malagasy: tsarò), Mom (Persian: mâmân, Malagasy: mama), and Dad (Persian: bâbâ, Malagasy: dada).")
    
    # This print statement is for the final answer block as requested by the prompt
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    find_similar_words()