def solve_problem():
    """
    This function solves the problem by calculating the sum of numbers in two given expression strings.

    Based on the analysis of the provided RASPy code, the main function 'v' checks a condition
    based on the function 'u'. For the given inputs, this condition is false.
    Therefore, the code proceeds to calculate the sum of the two numbers in the expression.
    This script performs that addition directly.
    """

    # First input expression
    expression1 = "734107+4295754"
    parts1 = expression1.split('+')
    num1_1 = int(parts1[0])
    num1_2 = int(parts1[1])
    result1 = num1_1 + num1_2

    # Second input expression
    expression2 = "5429141+142196"
    parts2 = expression2.split('+')
    num2_1 = int(parts2[0])
    num2_2 = int(parts2[1])
    result2 = num2_1 + num2_2

    # As per the instruction "Remember in the final code you still need to output
    # each number in the final equation!", we can show the equations.
    # However, the final format is specified as "output1;output2".
    # The printed output of this script will be in that final required format.
    
    # print(f"{num1_1} + {num1_2} = {result1}")
    # print(f"{num2_1} + {num2_2} = {result2}")

    # Print the final results in the required format "output1;output2".
    print(f"{result1};{result2}")

solve_problem()