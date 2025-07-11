def solve_riddle():
    """
    This function solves the riddle by combining clues from the image and the text.

    Clue 1: "X Y is pictured at an event on June 20, 2019."
    - The event is the 2019 NHL Awards.
    - The player, Nikita Kucherov, won three major awards, a feat known as a "hat trick".
    - This points to the word "hat".

    Clue 2: "Another X Y has almost ceased to be used due to the toxicity of mercury salts."
    - Mercury salts were used to make felt for hats.
    - This process was toxic, and the resulting object was a "felt hat".

    Solution: The pun connects the "hat" from "hat trick" to the "felt hat" from the mercury clue.
    """
    x = "Felt"
    y = "Hat"
    
    # The final answer is the combination of the two words.
    answer = f"{x} {y}"
    
    print(answer)

solve_riddle()