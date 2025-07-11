def solve_puzzle():
    """
    This function solves the video game trivia puzzle by finding the answers,
    extracting the second letter of each answer, and combining them to
    reveal the hidden word.
    """

    # The answers to the five trivia questions
    answer1 = "CRACKED"
    answer2 = "KONAMI"
    answer3 = "EARTHWORM"
    answer4 = "ADJUTANT"
    answer5 = "MOUSE AIM"

    answers = [answer1, answer2, answer3, answer4, answer5]

    # Extract the second letter from each answer
    hidden_letters = [ans[1] for ans in answers]

    # Combine the letters to form the hidden word
    hidden_word = "".join(hidden_letters)

    # Print the equation revealing the word
    print("The hidden word is formed by the second letters of each answer:")
    for i, answer in enumerate(answers):
        letter = answer[1]
        # Format the output to clearly show the source of each letter
        formatted_answer = f"({answer[0]})-{letter.upper()}-({answer[2:]})"
        print(f"({i+1}) Answer: {answer}, Second letter: {letter.upper()} -> Word part: {formatted_answer}")

    print("\nCombining the letters:")
    equation = " + ".join([f'"{letter}"' for letter in hidden_letters])
    print(f"{equation} = {hidden_word}")

solve_puzzle()

<<<ROADO>>>