def apply_t_operation(s):
    # Define the substrings and their corresponding insertions
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
    
    # Process the string from left to right
    while i <= len(s_list) - 4:
        for pattern, insert_char in operations:
            # Check if the current substring matches the pattern
            if s_list[i:i+4] == list(pattern):
                # Insert the character after the pattern
                s_list.insert(i+4, insert_char)
                # Move the index past the newly inserted character
                i += 5
                break
        else:
            # If no pattern matched, just move to the next character
            i += 1
    
    # Convert the list back to a string
    return ''.join(s_list)

# Given string
input_string = "BABCDDBCDECEABC"
# Apply the T operation
result = apply_t_operation(input_string)
print(result)