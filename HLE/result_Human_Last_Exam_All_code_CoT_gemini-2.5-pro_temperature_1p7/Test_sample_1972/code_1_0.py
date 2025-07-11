def find_borges_reference():
    """
    This function provides the answer to the user's literary question
    by storing the information in variables and printing it.
    """
    # Information about the reference
    novel_title = "El astillero (The Shipyard)"
    author_name = "Juan Carlos Onetti"
    source_of_quote = "his prologue to a 1968 edition of the novel"
    
    # The quotes from Borges
    quote_about_novel = "the intensity of a tiger and the variety that a chess duel can achieve"
    quote_about_author = "a continuator (and simplifier) of Faulkner"

    # Print the detailed answer
    print(f"The novel Jorge Luis Borges refers to is '{novel_title}' by {author_name}.")
    print("\nIn " + source_of_quote + ", Borges writes that the novel possesses \"" + quote_about_novel + "\".")
    print(f"He also describes its author, {author_name}, as {quote_about_author}.")

# Execute the function to display the information
find_borges_reference()