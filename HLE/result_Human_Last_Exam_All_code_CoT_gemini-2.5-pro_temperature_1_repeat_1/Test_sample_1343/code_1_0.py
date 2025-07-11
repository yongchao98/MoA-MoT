def solve_addition():
    """
    This function solves the problem by directly performing the addition operations
    described in the input strings. This is based on the analysis that the complex
    RASPy code is an obfuscated implementation of integer addition.
    """
    
    input1 = "734107+4295754"
    input2 = "5429141+142196"

    # Process the first input
    parts1 = input1.split('+')
    num1_1 = int(parts1[0])
    num1_2 = int(parts1[1])
    result1 = num1_1 + num1_2

    # Process the second input
    parts2 = input2.split('+')
    num2_1 = int(parts2[0])
    num2_2 = int(parts2[1])
    result2 = num2_1 + num2_2

    # The prompt also says to "output each number in the final equation".
    # This is slightly ambiguous, but we will interpret it as showing the work
    # before printing the final, machine-readable answer.
    # However, the top-level instruction says to only have one code block and
    # one response, so we will format the final output as requested.
    
    # Print the final result in the format "output1;output2"
    print(f"{result1};{result2}")

solve_addition()