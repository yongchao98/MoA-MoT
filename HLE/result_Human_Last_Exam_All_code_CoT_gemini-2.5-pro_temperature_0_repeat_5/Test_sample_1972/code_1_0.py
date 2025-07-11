def find_borges_reference():
    """
    This function identifies and prints the novel and author
    referenced by Jorge Luis Borges with the given quotes.
    """
    novel_title = "The Invention of Morel (La invención de Morel)"
    author_name = "Adolfo Bioy Casares"
    context = "the prologue to the 1940 edition of the novel"

    print("Jorge Luis Borges was referring to the novel:")
    print(f"'{novel_title}' by {author_name}.")
    print("\nHe made these comments in his famous " + context + ".")
    print("In that prologue, he praises the novel's plot, describing it as having 'the intensity of a tiger and the variety that a chess duel can achieve.'")
    print("He also positions Bioy Casares's clear, classical style as a continuation and simplification of the complexities found in authors like Faulkner.")

if __name__ == "__main__":
    find_borges_reference()