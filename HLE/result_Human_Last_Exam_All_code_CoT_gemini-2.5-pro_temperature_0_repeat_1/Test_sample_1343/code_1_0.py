def solve_addition():
    """
    This function solves the addition problems based on the analysis of the RASPy code.
    It calculates the sum for each input string and prints the results as a single
    semicolon-separated string.
    """
    # First input expression
    input1 = "734107+4295754"
    parts1 = input1.split('+')
    num1_1 = int(parts1[0])
    num1_2 = int(parts1[1])
    result1 = num1_1 + num1_2

    # Second input expression
    input2 = "5429141+142196"
    parts2 = input2.split('+')
    num2_1 = int(parts2[0])
    num2_2 = int(parts2[1])
    result2 = num2_1 + num2_2

    # The prompt asks to output each number in the final equation.
    # The following lines are for clarity, showing the individual calculations.
    # print(f"{num1_1} + {num1_2} = {result1}")
    # print(f"{num2_1} + {num2_2} = {result2}")

    # The final answer should be in the format "output1;output2".
    final_output = f"{result1};{result2}"
    print(final_output)

solve_addition()