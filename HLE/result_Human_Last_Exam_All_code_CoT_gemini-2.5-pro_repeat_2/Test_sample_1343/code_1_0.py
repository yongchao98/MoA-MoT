def solve_addition_puzzles():
    """
    This function solves two addition problems presented as strings,
    and prints the results in the specified format.
    """
    # First input string and calculation
    input1_str = "734107+4295754"
    parts1 = input1_str.split('+')
    num1_1 = int(parts1[0])
    num1_2 = int(parts1[1])
    result1 = num1_1 + num1_2

    # Second input string and calculation
    input2_str = "5429141+142196"
    parts2 = input2_str.split('+')
    num2_1 = int(parts2[0])
    num2_2 = int(parts2[1])
    result2 = num2_1 + num2_2

    # The prompt asks for the output of the two operations in the format "output1;output2"
    # and to output each number in the final equation. The code above calculates
    # the components of the equation (num1, num2, result), and the print statement
    # below formats the final output as requested.
    
    # Printing the first result
    print(result1, end=';')
    
    # Printing the second result
    print(result2)

solve_addition_puzzles()