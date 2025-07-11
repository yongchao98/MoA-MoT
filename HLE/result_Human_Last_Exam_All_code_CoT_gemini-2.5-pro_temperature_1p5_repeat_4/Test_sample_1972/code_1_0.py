def find_borges_reference():
    """
    This function identifies and prints the novel and author
    referenced by Jorge Luis Borges with the given quotes.
    """
    novel_title_original = "El astillero"
    novel_title_english = "The Shipyard"
    author = "Juan Carlos Onetti"
    source_of_quote = "Borges's prologue to the novel"

    print("The novel Jorge Luis Borges was referring to is:")
    print(f"- Title: {novel_title_original} ({novel_title_english})")
    print(f"- Author: {author}")
    print("\nBorges wrote these words in his prologue for the novel, where he also described Onetti's style in relation to William Faulkner's.")

if __name__ == "__main__":
    find_borges_reference()