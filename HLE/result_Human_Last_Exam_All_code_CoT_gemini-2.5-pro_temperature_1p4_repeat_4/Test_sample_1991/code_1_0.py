def solve_task():
    """
    This function generates a phrase according to the user's specifications
    and provides a syllable breakdown as a mathematical equation.
    """
    # The fourteen-syllable phrase about abnormal cell death causing cancer.
    # It is capitalized and its scansion matches "Supercalifragilisticexpialidocious".
    phrase = "APOPTOSIS FAILURE STARTS GRAVE TUMOROUS NECROSIS"

    # A dictionary mapping each word to its pre-verified syllable count.
    # This is a reliable method for demonstrating the syllable math.
    syllable_counts = {
        'APOPTOSIS': 4,
        'FAILURE': 2,
        'STARTS': 1,
        'GRAVE': 1,
        'TUMOROUS': 3,
        'NECROSIS': 3
    }

    words = phrase.split()
    
    # Retrieve the syllable count for each word from the dictionary.
    counts = [syllable_counts[word] for word in words]
    
    # Calculate the total syllables.
    total_syllables = sum(counts)

    # Format the list of counts into a string representing the equation.
    # This fulfills the requirement to "output each number in the final equation".
    equation_string = " + ".join(map(str, counts)) + f" = {total_syllables}"

    print("Syllable Equation for the Phrase:")
    print(equation_string)
    
    print("\nThe Final Phrase:")
    print(phrase)

solve_task()