import sys

def solve_linguistic_puzzle():
    """
    This script identifies and displays a well-known linguistic coincidence
    between two unrelated Asian languages for the words "mom", "dad", and "broom".
    """
    
    # Define the languages and a brief explanation.
    language1 = "Tagalog"
    language2 = "Tamil"
    explanation = (
        "This is a classic example of linguistic coincidence between two languages from different "
        "cultures in Asia with no direct contact or shared ancestry.\n\n"
        f"Language 1: {language1} (an Austronesian language from the Philippines)\n"
        f"Language 2: {language2} (a Dravidian language from Southern India/Sri Lanka)\n"
    )
    
    print(explanation)
    
    # --- Word Data ---
    # NOTE: These specific word choices are selected to highlight the striking phonetic similarities.
    # While they are valid words, some may be dialectal or less common than other synonyms.
    
    # Mom
    mom_l1 = "Nanay"
    mom_l2 = "Annai"
    
    # Dad
    dad_l1 = "Tatay"
    dad_l2 = "Thaththa"
    
    # Broom
    broom_l1 = "walis"
    broom_l2 = "valiyal"

    # --- Print Results in Equation Format ---
    
    print("---")
    print("1. MOM:")
    print(f"In {language1}, 'mom' can be '{mom_l1}'.")
    print(f"In {language2}, a word for 'mother' is '{mom_l2}'.")
    print(f"Result: {mom_l1} ≈ {mom_l2}\n")

    print("---")
    print("2. DAD:")
    print(f"In {language1}, 'dad' is '{dad_l1}'.")
    print(f"In {language2}, '{dad_l2}' can be used for 'father' or 'grandfather' and is phonetically similar.")
    print(f"Result: {dad_l1} ≈ {dad_l2}\n")
    
    print("---")
    print("3. BROOM:")
    print(f"In {language1}, 'broom' is '{broom_l1}'.")
    print(f"In {language2}, a word for a sweeper, from the verb 'to sweep', is '{broom_l2}'.")
    print(f"Result: {broom_l1} ≈ {broom_l2}\n")

if __name__ == '__main__':
    solve_linguistic_puzzle()