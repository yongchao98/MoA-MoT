import sys
import codecs

# Set stdout to utf-8 to ensure Chinese characters are printed correctly.
sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

def solve_riddle():
    """
    This function solves the Chinese character riddle and prints the explanation.
    """
    
    # The riddle describes a character by its constituent parts (radicals).

    # Part 1 of the riddle: "One horizontal stroke, another horizontal stroke, after another"
    # This describes the left part of the target character.
    # It refers to the three horizontal strokes in the '王' (wáng) radical.
    left_part_description = "'One horizontal stroke, another horizontal stroke, after another'"
    left_part_radical = "王"

    # Part 2 of the riddle: "one vertical stroke, another vertical stroke, after another; one vertical on the left, one vertical on the right..."
    # This describes the right part of the target character.
    # It accurately depicts the character '冊' (cè), which is composed of several vertical strokes representing a bound book of bamboo slips.
    right_part_description = "'one vertical stroke, another vertical stroke, after another...'"
    right_part_radical = "冊"

    # Combining the parts forms the final character.
    # Left part (王) + Right part (冊) = 珊 (shān)
    final_character = "珊"
    meaning = "coral"
    
    print("The riddle is solved by breaking it down into descriptions of the character's components:")
    print("-" * 20)
    
    print(f"1. {left_part_description}")
    print(f"   This describes the left-side radical: {left_part_radical}\n")
    
    print(f"2. {right_part_description}")
    print(f"   This describes the right-side component: {right_part_radical}\n")
    
    print(f"Combining these two parts [{left_part_radical} + {right_part_radical}] gives the final character.")
    print("-" * 20)
    
    print(f"The character is: {final_character}")
    print(f"Pinyin: shān")
    print(f"Meaning: {meaning}")

solve_riddle()