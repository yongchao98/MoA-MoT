def solve_wellington_mcq():
    """
    This function identifies the correct statements about the Duke of Wellington's career.
    Each statement is evaluated based on historical facts provided in the problem description.
    
    1. Correct: Indian logistical experience was crucial and adapted for the Peninsular War.
    2. Incorrect: His intelligence system was famously successful in the Peninsula (a European theater).
    3. Incorrect: The army did not standardize his structures in 1815; it largely reverted to older systems.
    4. Incorrect: Sandhurst was founded much earlier than 1829.
    5. Incorrect: Logistical lessons from India were highly relevant elsewhere.
    6. Correct: Integrating local forces was a key practice learned in India and used throughout the empire.
    7. Incorrect: Indian experience was the primary influence on his Peninsular logistics, not European theory.
    8. Correct: The "flying column" tactic was transferred from India to the Peninsula and later colonial wars.
    9. Incorrect: The 1813 Charter Act was driven by economic and political factors, not his military principles.
    10. Incorrect: His approach to civilian administration in India was a direct model for his actions in the Peninsula.
    
    The correct options are 1, 6, and 8.
    """
    
    correct_options = [1, 6, 8]
    
    # Sort the options in ascending order
    correct_options.sort()
    
    # Format the output as a comma-separated string
    # The problem asks to output each number in the final equation, which is interpreted as printing the final answer string.
    final_answer = ",".join(map(str, correct_options))
    
    print(final_answer)

solve_wellington_mcq()