def solve_riddle():
    """
    This script solves the riddle by identifying the two-word phrase 'X Y'.
    
    Clue 1: 'X Y' is pictured (represented) at the 2019 NHL Awards. In hockey, a 'hat trick'
              is a major achievement celebrated by fans throwing hats onto the ice.
    
    Clue 2: 'Another X Y' refers to the physical object, the 'Felt Hat', which was historically
              made using toxic mercury salts, a practice now ceased.
    """
    
    # The two words of the answer
    word_x = "Felt"
    word_y = "Hat"
    
    # Combine the words to form the final answer
    answer = f"{word_x} {word_y}"
    
    print(f"The answer to the riddle 'X Y' is: {answer}")

solve_riddle()