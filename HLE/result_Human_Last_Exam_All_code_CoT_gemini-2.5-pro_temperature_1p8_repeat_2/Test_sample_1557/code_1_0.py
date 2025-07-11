def solve_wellington_mcq():
    """
    This function identifies and prints the correct options for the historical question.
    Each statement is evaluated based on historical facts about the Duke of Wellington's career.
    
    1. Wellington's commissariat system from India was adapted for the Peninsula. (True)
    2. Intelligence system model for India but never used in Europe. (False - used effectively in Peninsula)
    3. Deccan force structure standardized across the army. (False - overstatement)
    4. Wellington established Sandhurst in 1829. (False - wrong date and context)
    5. Indian logistics had no impact elsewhere. (False - principles were transferable)
    6. Integrating local forces became standard colonial practice. (True)
    7. Peninsular reforms showed no Indian influence. (False - direct contradiction of facts)
    8. Flying columns adapted for Peninsula and Burma. (True)
    9. EIC reforms of 1813 influenced by Wellington's military principles. (False - drivers were political/economic)
    10. Indian civil administration approach had no bearing on the Peninsula. (False - direct continuity)
    
    The correct options are 1, 6, and 8.
    """
    
    # List of the numbers corresponding to the correct statements
    correct_options = [1, 6, 8]
    
    # Sort the numbers in ascending order (they are already sorted, but this is good practice)
    correct_options.sort()
    
    # Convert the numbers to strings to join them with a comma
    answer_string = ",".join(map(str, correct_options))
    
    # Print the final result
    print(answer_string)

solve_wellington_mcq()