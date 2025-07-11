def solve_riddle():
    """
    This function solves the Chinese character riddle by analyzing its parts.
    """
    # The riddle describes a character with multiple attributes.
    # Clue 1: "One horizontal stroke, another horizontal stroke, after another"
    # This suggests the character 三 (sān), which is what the answer looks like when rotated.
    clue_1_char = "三"

    # Clues 2, 3, and 4 point to the final answer.
    # Clue 2 & 4: "one vertical stroke, another vertical stroke, after another"
    # Clue 3: "one vertical on the left, one vertical on the right"
    # These clues describe the character 川 (chuān).
    final_answer_char = "川"
    
    # Explanation
    explanation = (
        "The riddle describes the character '川' (chuān), which means 'river'.\n\n"
        "Here is the breakdown:\n"
        "1. '...one vertical stroke, another vertical stroke, after another' (repeated) and 'one vertical on the left, one vertical on the right' perfectly describe the three vertical strokes of 川.\n"
        "2. 'One horizontal stroke, another horizontal stroke, after another' describes the character 三 (sān). This is a classic riddle technique, as 川 looks like 三 when rotated 90 degrees.\n\n"
        "Therefore, the character is:"
    )

    print(explanation)
    print(final_answer_char)

solve_riddle()