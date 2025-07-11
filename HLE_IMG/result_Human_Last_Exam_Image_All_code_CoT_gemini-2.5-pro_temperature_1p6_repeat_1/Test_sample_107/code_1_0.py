import sys
import io

# Helper function to emulate print for the final answer block
def new_print(*args, **kwargs):
    # If the string '<<<' is in the output, we'll know it's the final answer
    if '<<<' in args[0]:
        # This will be captured by the calling environment
        original_stdout.write(args[0] + '\n')
    else:
        # Regular prints go to the captured output
        captured_output.write(' '.join(map(str, args)) + '\n')

# This part is to make the explanation readable in the final output block
# It's a bit of a trick for the platform this code runs on.
original_stdout = sys.stdout
sys.stdout = new_print
captured_output = io.StringIO()


def solve_shogi_castle():
    """
    Identifies the Shogi castle from the image and options.
    """
    # Step 1: Analyze the visual information.
    # The image shows a defensive formation (a "castle") in Shogi.
    # Key pieces are the King (玉), two Silver Generals (銀), and two Gold Generals (金).
    # The King is in the center, flanked by the two Silvers on the same rank.
    # The two Golds are one rank behind and further to the sides of the Silvers.
    # This is a strong, symmetrical, central castle formation.
    
    # Step 2: Compare with known Shogi castles from the options.
    # - Mino (Q), Anaguma (F), Fortress/Yagura (C) are incorrect as they are side castles.
    # - Central House (O) is a general term for any castle with the King in the center.
    # - This specific robust central formation is famously associated with the legendary
    #   player Makoto Nakahara, who mastered central King strategies.
    #   Therefore, it is known as the "Nakahara" castle.

    # Step 3: Conclude the most accurate answer.
    # While it is a type of "Central House", its specific and famous name is "Nakahara".
    answer_option = 'H'
    answer_name = 'Nakahara'
    
    # Printing the explanation.
    print("Thinking Process:")
    print("1. The image displays a Shogi castle with the King (玉) in the central file.")
    print("2. The King is protected by two Silvers (銀) to its sides and two Golds (金) behind and to the sides.")
    print("3. This symmetrical and very strong central formation is known as the Nakahara castle (中原囲い), named after the famous player Makoto Nakahara who popularized it.")
    print(f"4. Comparing this with the options, option '{answer_option}' is '{answer_name}'.")

    # This part captures the text for the final explanation block
    explanation_text = captured_output.getvalue()
    sys.stdout = original_stdout
    print(explanation_text)
    
    # Final answer in the required format.
    print(f'<<<{answer_option}>>>')

solve_shogi_castle()