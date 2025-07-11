def solve_wellington_mcq():
    """
    This function identifies the correct statements regarding the Duke of Wellington's career.
    After analyzing all ten options based on historical records, the correct ones are:
    1: Wellington's mobile commissariat system from India was adapted for the Peninsular Campaign.
    6: Wellington's method of integrating local auxiliary forces became standard colonial practice.
    8: Wellington's system of deploying flying columns was adapted for the Peninsula and later colonial wars.
    """

    # The numbers of the correct statements
    correct_options = [1, 6, 8]
    
    # Sort the numbers in ascending order
    correct_options.sort()
    
    # Format the output as a comma-separated string
    # The problem statement says: "Remember in the final code you still need to output each number in the final equation!"
    # This is interpreted as clearly printing the final list of numbers.
    final_answer = ",".join(map(str, correct_options))
    
    print(final_answer)

solve_wellington_mcq()