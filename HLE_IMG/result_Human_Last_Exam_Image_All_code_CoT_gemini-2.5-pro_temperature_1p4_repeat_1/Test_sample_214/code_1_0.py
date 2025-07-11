def solve_poetry_form():
    """
    Analyzes syllable counts of given poetic lines to determine the poetic form.
    """
    # Define the line from the final poem in the sequence
    final_poem_line = "nacreous wavers"
    final_poem_words = final_poem_line.split()
    # Syllable counts for each word in the final poem's line
    final_poem_syllables = [3, 2] # na-cre-ous (3), wa-vers (2)

    # Define the line from the poem in the image
    image_line = "ghostly velum forms like a dance"
    # The words, treating "like a" as one unit due to poetic elision (synalepha)
    image_line_words = ["ghostly", "velum", "forms", "like a", "dance"]
    # Syllable counts for each word in the image line, applying elision to "like a"
    image_line_syllables = [2, 2, 1, 1, 1] # ghost-ly (2), ve-lum (2), forms (1), like a (1), dance (1)

    # --- Analysis of the first line ---
    total_syllables_final = sum(final_poem_syllables)
    print("Analysis of the line from the final poem: '{}'".format(final_poem_line))
    print("This line fits the 5-syllable requirement for a haiku's first or third line.")
    # Build and print the equation for the final poem's line
    equation_final = " + ".join(map(str, final_poem_syllables))
    print(f"Syllable equation: {equation_final} = {total_syllables_final}\n")

    # --- Analysis of the second line ---
    total_syllables_image = sum(image_line_syllables)
    print("Analysis of the line from the image: '{}'".format(image_line))
    print("Using poetic elision for 'like a' (1 syllable), this line fits the 7-syllable middle line of a haiku.")
    # Build and print the equation for the image line
    equation_image = " + ".join(map(str, image_line_syllables))
    print(f"Syllable equation: {equation_image} = {total_syllables_image}\n")
    
    print("Conclusion: The 5 and 7 syllable lines strongly indicate the poems follow the Haiku form.")

solve_poetry_form()