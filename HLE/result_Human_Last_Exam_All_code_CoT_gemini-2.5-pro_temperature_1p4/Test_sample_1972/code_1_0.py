def find_borges_reference():
    """
    Identifies and prints the novel and author Borges referred to
    with the given descriptions.
    """
    novel_title = "The Invention of Morel (La invención de Morel)"
    author_name = "Adolfo Bioy Casares"
    
    quote_1 = "the intensity of a tiger and the variety that a chess duel can achieve"
    description_2 = "a continuator (and simplifier) of Faulkner"

    print("Jorge Luis Borges was referring to the novel:")
    print(f"'{novel_title}' by {author_name}.")
    print("\nIn his prologue to the novel, he praised it for having...")
    print(f"...\"{quote_1}\"")
    print("\nHe also described its author as...")
    print(f"...\"{description_2}\".")

find_borges_reference()