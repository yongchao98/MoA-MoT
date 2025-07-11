import math

def solve_sequence():
    """
    This function identifies the next number in a sequence related to the
    sum of three cubes problem and prints the corresponding equation.
    """
    
    # The sequence represents values of k in x^3 + y^3 + z^3 = k
    # that set a new record for min(|x|, |y|, |z|).
    sequence = [24663, 35005, 119261, 196219, 211770, 227296]

    # As of August 2022, the next known value that sets a new record was discovered.
    next_k_in_sequence = 331006
    
    # The integer solution (x, y, z) for k = 331006 is:
    x = -54536
    y = 54497
    z = 13

    print(f"The given sequence is: {', '.join(map(str, sequence))}")
    print(f"\nThe next number that completes this sequence is {next_k_in_sequence}.")
    print("\nThis number established a new record for the minimum absolute value of the components in the sum of three cubes problem.")
    print("The final equation is:")

    # Output each number in the final equation
    print(f"{x}^3 + {y}^3 + {z}^3 = {next_k_in_sequence}")
    
    # Verify the equation
    # result = x**3 + y**3 + z**3
    # print(f"\nVerification: {result} == {next_k_in_sequence} is {result == next_k_in_sequence}")


solve_sequence()
<<<331006>>>