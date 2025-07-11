# This script identifies and prints the numbers of the correct statements
# from the historical analysis of the Duke of Wellington's career.

def solve_history_mcq():
    """
    Based on a detailed analysis of the Duke of Wellington's career, this function
    identifies the correct statements regarding the impact of his innovations.

    - Statement 1: Correct. Indian logistics were adapted for the Peninsular War.
    - Statement 6: Correct. The use of local auxiliary forces became standard practice.
    - Statement 8: Correct. The "flying column" tactic was adapted for other campaigns.
    
    The other statements are historically inaccurate.
    """
    
    # A list containing the numbers of the correct options.
    correct_options = [1, 6, 8]
    
    # The numbers should be in sorted order.
    correct_options.sort()
    
    # The final output is a string of the numbers separated by a comma.
    final_answer = ",".join(map(str, correct_options))
    
    print(final_answer)

solve_history_mcq()