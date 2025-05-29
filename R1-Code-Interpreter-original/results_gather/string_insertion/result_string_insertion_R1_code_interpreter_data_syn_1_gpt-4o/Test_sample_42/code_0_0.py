def apply_t_operation(s):
    # Define the patterns and their corresponding insertions
    patterns = [
        ("ABCD", "A"),
        ("BCDE", "B"),
        ("CDEA", "C"),
        ("DEAB", "D"),
        ("EABC", "E")
    ]
    
    # Convert the string to a list for easier manipulation
    s_list = list(s)
    i = 0
    
    # Iterate over the string
    while i <= len(s_list) - 4:
        # Check each pattern
        for pattern, insert_char in patterns:
            # If the pattern matches, insert the character
            if s_list[i:i+4] == list(pattern):
                s_list.insert(i+4, insert_char)
                i += 4  # Move past the inserted character
                break
        else:
            i += 1  # Move to the next character if no pattern matched
    
    # Convert the list back to a string
    return ''.join(s_list)

# Given string
initial_string = "EABCCEABCACDEA"

# Apply the T operation
result_string = apply_t_operation(initial_string)

# Print the result
print(result_string)