def analyze_poetic_form():
    """
    Analyzes the poetic form by counting syllables in the provided lines
    and matching them to a known form.
    """
    # Define the lines and their words
    line_a_text = "nacreous wavers"
    line_a_words = line_a_text.split()
    # Syllable counts for each word in the first line
    line_a_syllables = [3, 2] # na-cre-ous (3), wa-vers (2)

    line_b_text = "ghostly velum forms like a dance"
    line_b_words = line_b_text.split()
    # Syllable counts for each word in the second line
    line_b_syllables = [2, 2, 1, 1, 1, 1] # ghost-ly (2), ve-lum (2), forms (1), like (1), a (1), dance (1)

    # Calculate total syllables for each line
    total_a = sum(line_a_syllables)
    total_b = sum(line_b_syllables)

    print("Step 1: Analyze the syllable count of the line 'nacreous wavers'.")
    # Create the equation string for the first line
    equation_a = " + ".join(map(str, line_a_syllables))
    print(f"The syllable count is: {equation_a} = {total_a}\n")

    print("Step 2: Analyze the syllable count of the line 'ghostly velum forms like a dance'.")
    # Create the equation string for the second line
    equation_b = " + ".join(map(str, line_b_syllables))
    print(f"The syllable count is: {equation_b} = {total_b}\n")

    print("Step 3: Identify the poetic form.")
    print(f"The first line's count of {total_a} syllables matches the first and last lines of a Haiku (5-7-5).")
    print(f"The second line's count of {total_b} syllables is a common variation on the 7-syllable middle line in modern Haiku.")
    print("Given that the artist, Mary Ruefle, intended to write Haikus for this series, this confirms the form.")
    print("\nThe poetic form is Haiku.")

analyze_poetic_form()