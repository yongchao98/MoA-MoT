def solve_addition():
    """
    This function solves the addition problems presented in the two input strings
    and prints the results in the specified format.
    """
    # First input string
    input1 = "734107+4295754"
    num1_str, num2_str = input1.split('+')
    num1 = int(num1_str)
    num2 = int(num2_str)
    result1 = num1 + num2

    # Second input string
    input2 = "5429141+142196"
    num3_str, num4_str = input2.split('+')
    num3 = int(num3_str)
    num4 = int(num4_str)
    result2 = num3 + num4

    # Print the outputs in the "output1;output2" format
    print(f"{result1};{result2}")

solve_addition()