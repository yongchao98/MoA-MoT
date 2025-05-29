def apply_t_operation(s):
    # Define the operations as a list of tuples (substring, character to insert)
    operations = [
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
    while i < len(s_list) - 3:
        # Check each operation
        for sub, char in operations:
            # If the substring is found, insert the character
            if s_list[i:i+4] == list(sub):
                s_list.insert(i+4, char)
                i += 4  # Move past the inserted character
                break  # Break to start checking from the new position
        else:
            i += 1  # Move to the next character if no operation was performed
    
    # Convert the list back to a string
    return ''.join(s_list)

# Given string
initial_string = "EAABCDCABCDADEAB"
# Apply the T operation
result_string = apply_t_operation(initial_string)
print(result_string)