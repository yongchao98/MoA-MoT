def apply_t_operation(s):
    # Define the rules for the T operation
    rules = [
        ("ABCD", "A"),
        ("BCDE", "B"),
        ("CDEA", "C"),
        ("DEAB", "D"),
        ("EABC", "E")
    ]
    
    # Convert the string to a list for easier manipulation
    s = list(s)
    i = 0
    
    # Iterate over the string
    while i < len(s) - 3:
        # Check each rule
        for pattern, insert_char in rules:
            # If the pattern matches, insert the character
            if s[i:i+4] == list(pattern):
                s.insert(i+4, insert_char)
                i += 4  # Move past the inserted character
                break  # Restart checking from the current position
        i += 1
    
    # Convert the list back to a string
    return ''.join(s)

# Given string
initial_string = "CDEAEBCDECCDEAED"

# Apply the T operation
result_string = apply_t_operation(initial_string)

# Print the result
print(result_string)