def solve_trivia_puzzle():
    """
    Solves the trivia puzzle by identifying the answer to each question
    and combining the first letters to form a hidden word.
    """

    # Answers to the trivia questions
    answer1 = "Slumdog Millionaire"
    answer2 = "Hell"
    answer3 = "ALF"
    answer4 = "Winnie-the-Pooh"

    # Extract the first letter of each answer
    first_letters = [answer1[0], answer2[0], answer3[0], answer4[0]]

    # Combine the letters to form the hidden word
    hidden_word = "".join(first_letters)

    # Print the explanations for each answer
    print("(1) The sequence of multiple-choice answers (A, A, D, A, C, A, D, D, A) corresponds to Jamal's path to victory in the Oscar-winning film 'Slumdog Millionaire'.")
    print(f"The first letter is {first_letters[0]}.\n")

    print("(2) In medieval Christian diagrams of the cosmos, Hell is traditionally depicted at the bottom. The clue about Stalin refers to screen direction for 'bad guys', fitting the concept of Hell.")
    print(f"The first letter is {first_letters[1]}.\n")

    print("(3) The cat-themed dishes 'Cheshire Salad' and 'Siamese Hotpot' point to the 1980s TV series 'ALF', whose titular character frequently tried to eat cats.")
    print(f"The first letter is {first_letters[2]}.\n")

    print("(4) The character 'Winnie-the-Pooh' has been banned in China due to memes comparing his appearance to leader Xi Jinping. The puzzle uses 'X' to refer to both the character and his media sphere.")
    print(f"The first letter is {first_letters[3]}.\n")

    # Print the final equation and the hidden word
    print("The hidden word is formed by the first letter of each answer:")
    final_equation = f"{first_letters[0]} + {first_letters[1]} + {first_letters[2]} + {first_letters[3]} = {hidden_word}"
    print(final_equation)

solve_trivia_puzzle()