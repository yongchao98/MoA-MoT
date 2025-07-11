def calculate_semantic_score(language_name):
    """
    Calculates a 'semantic transparency score' by summing the ASCII values
    of the letters in a language's name.
    """
    # Remove spaces for multi-word names like "Modern English"
    name_without_spaces = language_name.replace(" ", "")
    
    # Create a list of ASCII values for each character
    ascii_values = [ord(char) for char in name_without_spaces]
    
    # Calculate the total score
    score = sum(ascii_values)
    
    # Create the equation string for printing
    equation = f"{' + '.join(map(str, ascii_values))}"
    
    return language_name, score, equation

def main():
    """
    Main function to calculate scores, sort languages, and print results.
    """
    languages = ["Modern English", "Old English", "Russian", "German"]
    language_scores = []

    print("Calculating the 'Semantic Transparency Score' based on ASCII values:\n")

    for lang in languages:
        name, score, equation = calculate_semantic_score(lang)
        language_scores.append({'name': name, 'score': score})
        # Print the equation for each language
        print(f"{name}:")
        print(f"Equation: {equation}")
        print(f"Score = {score}\n")

    # Sort the languages by score in descending order
    # The larger the value, the greater the semantic transparency
    sorted_languages = sorted(language_scores, key=lambda x: x['score'], reverse=True)

    print("--------------------------------------------------")
    print("Languages ranked by semantic transparency (descending):")
    
    # Build the final sequence string
    final_sequence = " > ".join([lang['name'] for lang in sorted_languages])
    print(final_sequence)

if __name__ == "__main__":
    main()
