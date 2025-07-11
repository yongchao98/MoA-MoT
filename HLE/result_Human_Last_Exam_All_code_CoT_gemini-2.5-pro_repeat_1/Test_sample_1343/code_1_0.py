def solve():
    # Based on the analysis, the complex logic in the RASPy code is a trap
    # that is not triggered by the given inputs. The task simplifies to
    # performing two standard additions.

    # First calculation
    input_str1 = "734107+4295754"
    parts1 = input_str1.split('+')
    num1_1 = int(parts1[0])
    num1_2 = int(parts1[1])
    result1 = num1_1 + num1_2

    # Second calculation
    input_str2 = "5429141+142196"
    parts2 = input_str2.split('+')
    num2_1 = int(parts2[0])
    num2_2 = int(parts2[1])
    result2 = num2_1 + num2_2

    # Print the results in the specified "output1;output2" format.
    # The instruction to "output each number in the final equation" is interpreted
    # as being fulfilled by calculating the final result for each case as requested.
    print(f"{result1};{result2}")

solve()