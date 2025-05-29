def perform_t_operation(s):
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
    
    # Iterate over the string
    while i < len(s_list) - 3:
        # Check each operation
        for pattern, insert_char in operations:
            # If the pattern is found, insert the character
            if s_list[i:i+4] == list(pattern):
                s_list.insert(i+4, insert_char)
                # Move the index forward to avoid re-checking the same pattern
                i += 4
                break
        else:
            # If no pattern matched, move to the next character
            i += 1
    
    # Convert the list back to a string
    result = ''.join(s_list)
    return f"<<<{result}>>>"

# Given string
input_string = "ABCDDCDEAEDEABD"
# Perform the operation and print the result
print(perform_t_operation(input_string))