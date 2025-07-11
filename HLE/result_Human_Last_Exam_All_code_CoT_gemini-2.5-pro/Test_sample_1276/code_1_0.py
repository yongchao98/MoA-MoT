def solve_puzzle():
    """
    Solves the 'Fun Facts From The Zoo' puzzle.

    The puzzle's solution is derived by interpreting the flavor text "Why was the earliest
    missive so similar to 'A'?" as a clue to the extraction method. "Earliest missive"
    points to taking the first letter of the answers to the riddles. "Similar to 'A'"
    is a pun on the article "a", suggesting the sound of a letter, like "a 'W'".
    For a three-word phrase, we look for three such clues, leading to "WWW" or
    "WORLD WIDE WEB". The many other clues are red herrings.

    Here are the three key clues and their solutions:
    """

    print("Step 1: Identifying the first 'W'")
    print("Clue: Why was the baby’s crying violating the terms of a 1946 international convention? (5)")
    print("Solution: The convention is for regulating 'whaling'. The baby's crying is 'wailing'. The answer is WAILS.")
    first_word_clue = "WAILS"
    print(f"The first letter is {first_word_clue[0]}.\n")

    print("Step 2: Identifying the second 'W'")
    print("Clue: What do you call a snake that cleans windows? (5)")
    print("Solution: This is a pun on 'viper' and 'wiper'. The answer is WIPER.")
    second_word_clue = "WIPER"
    print(f"The second letter is {second_word_clue[0]}.\n")

    print("Step 3: Identifying the third 'W'")
    print("Clue: Why did the mercury-laden, overfished fish’s string instrument sound so bad? (4)")
    print("Solution: The fish is a 'tuna', a pun on 'tune'. The instrument sounded bad because the sound was WEAK.")
    third_word_clue = "WEAK"
    print(f"The third letter is {third_word_clue[0]}.\n")

    print("Step 4: Assembling the final answer")
    print("The three 'W's form the acronym 'WWW'.")
    final_phrase = ["WORLD", "WIDE", "WEB"]
    print(f"The final phrase is: {final_phrase[0]} {final_phrase[1]} {final_phrase[2]}")

solve_puzzle()
<<<WORLD WIDE WEB>>>