def apply_T_operation(s):
    # Define the substrings and their corresponding insertions
    operations = [
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
    while i <= len(s) - 4:
        # Check each operation
        for pattern, insert_char in operations:
            # If the pattern matches, insert the character
            if s[i:i+4] == list(pattern):
                s.insert(i+4, insert_char)
                i += 4  # Move past the inserted character
                break
        else:
            i += 1  # Move to the next character if no pattern matched
    
    # Convert the list back to a string
    return ''.join(s)

# Given string
input_string = "BBEABCCABCDADEABB"
# Apply the T operation
result = apply_T_operation(input_string)
print(result)