def solve_kalabari_attire_question():
    """
    This function identifies and prints the correct answer to the multiple-choice question
    about the Kalabari 'etibo' shirt.
    """
    # The provided answer choices are stored in a dictionary.
    answer_choices = {
        'A': 'No special category',
        'B': 'Amanyanabo (sovereign)',
        'C': 'Alabo (chiefs)',
        'D': 'Asawo (young men)',
        'E': 'Opu asawo (gentlemen)'
    }

    # Based on cultural knowledge, the 'etibo' is a status symbol associated
    # with respected men of high standing, often referred to as gentlemen.
    # While chiefs (Alabo) would wear it, 'Opu asawo' is the broader and more
    # accurate category representing the class of men for whom this is characteristic attire.
    correct_answer_key = 'E'

    # Get the full description of the correct answer.
    correct_answer_description = answer_choices[correct_answer_key]

    print("Question: To which category of Kalabari men is the etibo dress associated?")
    print("-" * 20)
    print("Finding the correct association...")
    print(f"The 'etibo' is most closely associated with the Kalabari men of status and style known as 'Opu asawo'.")
    print(f"The correct choice is therefore '{correct_answer_key}'.")
    print("-" * 20)
    print(f"Final Answer: {correct_answer_key}. {correct_answer_description}")

# Run the function to display the answer.
solve_kalabari_attire_question()