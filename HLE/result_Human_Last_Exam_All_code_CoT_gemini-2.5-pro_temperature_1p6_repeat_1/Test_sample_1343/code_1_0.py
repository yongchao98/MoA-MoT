def solve_RASPy_logic(sop: str) -> str:
    """
    This function implements the distilled logic from the provided RASPy code.
    
    The core logic depends on a "trigger" condition: counting the occurrences
    of the digit '7' in three specific slices of the input string:
    - The first 3 characters (sop[0:3])
    - The characters at indices 3 through 5 (sop[3:6])
    - The last 3 characters (sop[-3:])

    If the total count of '7's is exactly 1, a special operation occurs.
    Otherwise, the function treats the input as a mathematical addition
    of two numbers separated by a '+' and returns their sum.
    """
    
    # This block corresponds to the logic in function u(sop)
    s1 = sop[0:3]
    s2 = sop[3:6]
    s3 = sop[-3:]
    
    count_of_sevens = s1.count('7') + s2.count('7') + s3.count('7')
    
    # This is the main branching logic from function v(sop)
    if count_of_sevens == 1:
        # This branch is not triggered by the given inputs.
        # The original code would produce a "get pwned" message.
        return "pwned"
    else:
        # This branch corresponds to the addition logic in v(sop).
        # We split the string by '+' and perform the addition.
        try:
            parts = sop.split('+')
            num1 = int(parts[0])
            num2 = int(parts[1])
            result = num1 + num2
            return str(result)
        except (ValueError, IndexError):
            return "Error: Invalid input format"

# The inputs provided in the problem statement
input_1 = "734107+4295754"
input_2 = "5429141+142196"

# Calculate the results for both inputs
output_1 = solve_RASPy_logic(input_1)
output_2 = solve_RASPy_logic(input_2)

# Print the final answer in the required "output1;output2" format
print(f"{output_1};{output_2}")