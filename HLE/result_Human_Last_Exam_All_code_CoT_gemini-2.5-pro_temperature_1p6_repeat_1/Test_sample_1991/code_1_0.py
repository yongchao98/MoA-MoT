def generate_cancer_phrase():
    """
    This function generates and prints a fourteen-syllable phrase
    about abnormal cell death, matching the scansion of a famous
    Disney song title. It also provides the syllable count equation.
    """

    # Define the words and their corresponding syllable counts
    words = ["When", "Subthalamic", "Cells", "Do", "Not", "Go", "Through", "Apoptosis"]
    syllables = [1, 4, 1, 1, 1, 1, 1, 4]

    # Capitalize each word (already done, but good practice for formatting)
    capitalized_words = [word.capitalize() for word in words]

    # Create the final phrase
    final_phrase = " ".join(capitalized_words)

    # Create the syllable count equation string
    equation_parts = [str(s) for s in syllables]
    total_syllables = sum(syllables)
    equation_string = f"Syllable Count: {' + '.join(equation_parts)} = {total_syllables}"

    # Print the results
    print(final_phrase)
    print(equation_string)

# Execute the function to display the output
generate_cancer_phrase()
