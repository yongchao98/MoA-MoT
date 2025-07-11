def solve_riddle():
    """
    Solves the four trivia questions and reveals the hidden word.
    """
    # Answers to the trivia questions
    answer1 = "Slumdog Millionaire"
    answer2 = "Orchestra"
    answer3 = "Red Dwarf"
    answer4 = "Topolino"

    # Extract the first letter of each answer
    letter1 = answer1[0]
    letter2 = answer2[0]
    letter3 = answer3[0]
    letter4 = answer4[0]

    # Combine letters to form the hidden word
    hidden_word = letter1 + letter2 + letter3 + letter4

    # Print the explanation for each letter
    print(f"(1) The sequence of answers A, A, D, A, C, A, D, D, A belongs to Jamal in '{answer1}'. The first letter is {letter1}.")
    print(f"(2) The '{answer2}' is at the bottom (in the pit) and its main section (violins) is on the left. The first letter is {letter2}.")
    print(f"(3) 'Cheshire Salad' and 'Siamese Hotpot' are cat-themed joke dishes from the sci-fi series '{answer3}'. The first letter is {letter3}.")
    print(f"(4) The puppet show '{answer4}' was banned in Argentina as its name resembled dictator Galtieri's nickname. The first letter is {letter4}.")
    print("-" * 20)
    print(f"The hidden word is: {hidden_word}")

solve_riddle()
<<<SORT>>>