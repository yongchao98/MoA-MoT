def solve_wellington_mcq():
    """
    This function identifies the correct statements about the Duke of Wellington's
    career based on historical analysis and prints the corresponding numbers in sorted order.

    - Statement 1: Correct. His Indian logistical system was adapted for the Peninsular War.
    - Statement 6: Correct. His methods for integrating local forces became standard colonial practice.
    - Statement 8: Correct. The "flying column" tactic from India was used in the Peninsula and later in Burma.
    - All other statements contain historical inaccuracies or unsubstantiated claims.
    """

    # The numbers of the statements identified as historically accurate.
    correct_options = [1, 6, 8]

    # Sort the list of numbers in ascending order.
    correct_options.sort()

    # The prompt requests the final answer to be the numbers separated by a comma.
    # We use the print function with a separator argument to format the output.
    # The '*' operator unpacks the list elements as individual arguments to the print function.
    print("The correct option numbers are:")
    print(*correct_options, sep=",")

solve_wellington_mcq()
#<<<1,6,8>>>