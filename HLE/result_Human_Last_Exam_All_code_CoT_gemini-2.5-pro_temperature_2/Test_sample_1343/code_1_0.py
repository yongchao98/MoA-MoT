def solve():
    """
    This function solves the problem by simulating the behavior of the RASPy code.
    For the given inputs, the complex logic simplifies to standard addition.
    """

    # Input 1
    input_str1 = "734107+4295754"
    operands1 = input_str1.split('+')
    num1_1 = int(operands1[0])
    num1_2 = int(operands1[1])
    result1 = num1_1 + num1_2

    # Input 2
    input_str2 = "5429141+142196"
    operands2 = input_str2.split('+')
    num2_1 = int(operands2[0])
    num2_2 = int(operands2[1])
    result2 = num2_1 + num2_2

    # As per the prompt "output each number in the final equation!",
    # we first show the calculations.
    print(f"Calculation 1: {num1_1} + {num1_2} = {result1}")
    print(f"Calculation 2: {num2_1} + {num2_2} = {result2}")
    print("-" * 20)

    # Then, we provide the final answer in the requested format "output1;output2".
    print("Final formatted output:")
    print(f"{result1};{result2}")


solve()