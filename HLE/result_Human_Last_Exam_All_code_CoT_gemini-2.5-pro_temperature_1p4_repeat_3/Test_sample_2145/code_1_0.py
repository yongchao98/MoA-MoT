def solve_trivia_puzzle():
    """
    This function solves a multi-part trivia question to find a hidden word.
    It prints the answer to each question, extracts the first letter,
    and constructs the final word.
    """

    # The answers to the four trivia questions
    answers = {
        1: "Amadeus",
        2: "Hell",
        3: "ALF",
        4: "Animal Farm"
    }

    # Brief explanations for each answer
    explanations = {
        1: "The sequence A, A, D, A C, A, D, D, A represents musical syllables (La, La, Re, La, Do, etc.) from a Mozart piece heavily featured in the Oscar-winning film 'Amadeus'.",
        2: "In filmmaking, antagonists or evil forces are often staged on the left. In historical and religious diagrams (e.g., Dante's Inferno), 'Hell' is depicted at the very bottom.",
        3: "The series is 'ALF'. The title character, an alien, famously wanted to eat cats. 'Cheshire Salad' and 'Siamese Hotpot' are humorous fan-invented dish names referencing this.",
        4: "The allegorical film 'Animal Farm' was widely banned. Its main characters, the pigs, are direct representations of totalitarian leaders, and the allegory applies to many dictators who rose through coups."
    }

    hidden_word = ""
    
    print("Let's solve the puzzle step by step:\n")

    # Process each answer
    for i in sorted(answers.keys()):
        answer = answers[i]
        first_letter = answer[0]
        hidden_word += first_letter
        
        print(f"Question ({i}):")
        print(f"The answer is '{answer}'.")
        print(f"Explanation: {explanations[i]}")
        print(f"The first letter is '{first_letter}'.\n")

    print("---")
    print("The hidden word is formed by concatenating the first letters of the answers:")
    
    # Show the formation of the word
    letters = " + ".join([f"'{ans[0]}'" for ans in answers.values()])
    print(f"{letters} = {hidden_word}")

solve_trivia_puzzle()