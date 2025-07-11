def solve_hidden_word():
    """
    This script solves the video game trivia puzzle to reveal a hidden word.
    It determines the answer to each trivia question, extracts the second letter,
    and combines them to form the final word.
    """
    
    # Answers to the five trivia questions
    answer1 = "OCCUPY"
    answer2 = "KONAMI"
    answer3 = "DESCENT"
    answer4 = "INSANUS"
    answer5 = "GAMING HAND"
    
    answers = [answer1, answer2, answer3, answer4, answer5]
    
    # Extract the second letter from each answer
    second_letters = [ans[1] for ans in answers]
    
    # The hidden word formed by the second letters
    hidden_word = "".join(second_letters)
    
    # Print the equation that reveals the hidden word
    print("The answers to the trivia questions are:")
    print(f"(1) {answer1}")
    print(f"(2) {answer2}")
    print(f"(3) {answer3}")
    print(f"(4) {answer4}")
    print(f"(5) {answer5}")
    print("\nCombining the second letter of each answer reveals the hidden word:")
    
    # Print the equation with each letter
    equation = " + ".join(second_letters)
    print(f"{equation} = {hidden_word}")

solve_hidden_word()
<<<COENA>>>