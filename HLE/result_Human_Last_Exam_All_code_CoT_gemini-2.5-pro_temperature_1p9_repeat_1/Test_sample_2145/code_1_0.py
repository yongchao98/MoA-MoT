def solve_trivia_puzzle():
    """
    Solves a multi-part trivia puzzle to reveal a hidden word.
    The final answer is formed by the first letter of each trivia answer.
    """

    # --- Question 1 ---
    answer1 = "Slumdog Millionaire"
    letter1 = answer1[0]
    print("(1) The Oscar-winning film where the main character's 'path to success' is shown through a sequence of game show answers is '{}'. The first letter is {}.".format(answer1, letter1))

    # --- Question 2 ---
    answer2 = "Devils"
    letter2 = answer2[0]
    print("(2) In film, 'enemies' or '{}' were placed on the left by Stalin's order. In classic iconography, they are depicted at the bottom (in Hell). The first letter is {}.".format(answer2, letter2))

    # --- Question 3 ---
    answer3 = "House"
    letter3 = answer3[0]
    print("(3) The TV series is '{}'. The clue 'Cheshire' points to Alice in Wonderland. Adding a 'missing letter' to 'Alice' creates 'Malice', hinting at the main character's personality. The first letter is {}.".format(answer3, letter3))

    # --- Question 4 ---
    answer4 = "Chipollino"
    letter4 = answer4[0]
    print("(4) The animation with characters resembling coup leaders is '{}'. The story is about an uprising of vegetables against fruit aristocrats. The first letter is {}.".format(answer4, letter4))

    # --- Final Hidden Word ---
    hidden_word = letter1 + letter2 + letter3 + letter4
    print("\nCombining the first letters gives the hidden word:")
    print("{} + {} + {} + {} = {}".format(letter1, letter2, letter3, letter4, hidden_word))

solve_trivia_puzzle()

# The final answer in the requested format
final_answer = "SDHC"
print("\n<<<{}>>>".format(final_answer))