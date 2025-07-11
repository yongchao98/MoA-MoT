import sys

def simulate_v(sop: str) -> str:
    """
    Simulates the behavior of the provided RASPy function `v` on a given input string `sop`.
    
    The function's behavior branches based on a condition checked by function `u`,
    which looks for a specific pattern of '7's in the input string.
    """
    
    # The logic of function `u` determines which path to take.
    # `u` checks for a pattern of '7's in three specific parts of the input string:
    # `c1`: Number of '7's in the first 3 characters (`sop[0:3]`)
    # `c2`: Number of '7's in the characters from index 3 to 5 (`sop[3:6]`)
    # `c3`: Number of '7's in the last 3 characters (`sop[-3:]`)
    c1 = sop[0:3].count('7')
    c2 = sop[3:6].count('7')
    c3 = sop[-3:].count('7')
    
    # After a complex sequence of operations, the function `v` checks if the result
    # of `u` is equal to 1. This condition is met if and only if:
    # `c1` is 1, AND `c2` is >= 1, AND `c3` is >= 1.
    is_pwned_trigger = (c1 == 1 and c2 >= 1 and c3 >= 1)

    if is_pwned_trigger:
        # Path 1: The "pwned" message is generated.
        # The output string must have the same length as the input `sop`.
        length = len(sop)
        base_msg = "get pwned"
        pwned_list = list(base_msg)
        
        # The rest of the string is filled with '!'s, with a special case for the
        # last character, which becomes '1' if its index is greater than 10.
        last_char_index = length - 1
        
        # Fill the remaining characters of the list
        for i in range(len(base_msg), length):
            if i == last_char_index and i > 10:
                pwned_list.append('1')
            else:
                pwned_list.append('!')
                
        return "".join(pwned_list)
        
    else:
        # Path 2: Standard addition is performed.
        # The complex RASPy functions for addition (`a`, `m`, `f`, `n`) result
        # in a standard sum of the two numbers.
        try:
            num1_str, num2_str = sop.split('+')
            num1 = int(num1_str)
            num2 = int(num2_str)
            result = num1 + num2
            
            # The prompt asks to omit leading zeros, which int-to-string conversion does automatically.
            return str(result)
        except (ValueError, IndexError):
            # This case shouldn't be reached with the given valid inputs.
            return "Error: Invalid input format for addition."

# The two inputs from the prompt
input1 = "734107+4295754"
input2 = "5429141+142196"

# Calculate the output for each input by simulating the RASPy code
output1 = simulate_v(input1)
output2 = simulate_v(input2)

# Print the final result in the specified "output1;output2" format
print(f"{output1};{output2}")