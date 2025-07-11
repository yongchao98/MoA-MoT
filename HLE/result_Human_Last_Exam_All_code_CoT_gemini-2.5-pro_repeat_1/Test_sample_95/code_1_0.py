def solve_riddle():
    """
    This function decodes the answer to the riddle using ASCII values.
    The riddle's clues point to The Alps:
    1. Visible from Milan but not from smoggy northern European cities.
    2. Requires 'imagination' to see them as more than just mountains (e.g., castles in the sky).
    """
    print("Decoding the answer from a numerical 'equation'...")
    print("-" * 40)

    # The answer is "The Alps". We'll use the ASCII numbers for each character.
    # First word: "The"
    word1_nums = [84, 104, 101]
    word1_chars = [chr(n) for n in word1_nums]
    word1 = "".join(word1_chars)

    # Second word: "Alps"
    word2_nums = [65, 108, 112, 115]
    word2_chars = [chr(n) for n in word2_nums]
    word2 = "".join(word2_chars)

    # Print the equation and numbers used for the first word
    print(f"Equation Part 1 uses numbers: {', '.join(map(str, word1_nums))}")
    print(f"Result -> chr({word1_nums[0]}) + chr({word1_nums[1]}) + chr({word1_nums[2]}) = '{word1}'")
    print()

    # Print the equation and numbers used for the second word
    print(f"Equation Part 2 uses numbers: {', '.join(map(str, word2_nums))}")
    print(f"Result -> chr({word2_nums[0]}) + chr({word2_nums[1]}) + chr({word2_nums[2]}) + chr({word2_nums[3]}) = '{word2}'")
    print("-" * 40)

    # Print the final combined answer
    final_answer = f"{word1} {word2}"
    print(f"The final answer to the riddle is: {final_answer}")

solve_riddle()