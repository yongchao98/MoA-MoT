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
    
    # Iterate over the string
    while i < len(s_list) - 3:
        for pattern, insert_char in operations:
            # Check if the current slice matches the pattern
            if s_list[i:i+4] == list(pattern):
                # Insert the character after the pattern
                s_list.insert(i+4, insert_char)
                # Move the index past the newly inserted character
                i += 4
                break
        i += 1
    
    # Convert the list back to a string
    result = ''.join(s_list)
    return f"<<<{result}>>>"

# Given string
input_string = "ABCDCCDEADCDEAEB"
# Apply the T operation
result = apply_t_operation(input_string)
# Print the result
print(result)