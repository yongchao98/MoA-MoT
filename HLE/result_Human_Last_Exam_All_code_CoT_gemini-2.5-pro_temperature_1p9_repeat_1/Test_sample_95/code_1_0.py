def solve_riddle():
    """
    This function decodes the riddle by connecting its historical and geographical clues
    to a famous 19th-century astronomical debate.
    """
    # The key year when Giovanni Schiaparelli made his observations from Milan.
    year_of_observation = 1877
    
    # The celestial body at the center of the debate.
    celestial_body = "Mars"
    
    # The name of the phenomenon, whose English translation caused the controversy.
    phenomenon = "Canals"
    
    # The riddle requires an "equation". We will construct a symbolic one
    # to represent the solution. The instruction is to output each number, so we will print it.
    print("Constructing the solution from the riddle's clues:")
    print("Key Number (Year of Observation):")
    print(year_of_observation)
    print("\nLeads to the identity of THEM, which are the...")
    
    # Combine the elements for the final answer.
    final_answer_part_1 = phenomenon
    final_answer_part_2 = "on"
    final_answer_part_3 = celestial_body
    
    print(f"{final_answer_part_1} {final_answer_part_2} {final_answer_part_3}")

solve_riddle()