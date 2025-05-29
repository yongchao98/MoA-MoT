def apply_t_operation(s):
    # Define the substrings and their corresponding insertions
    operations = [
        ("ABCD", "A"),
        ("BCDE", "B"),
        ("CDEA", "C"),
        ("DEAB", "D"),
        ("EABC", "E")
    ]
    
    # Initialize the result as a list of characters for easy insertion
    result = list(s)
    i = 0
    
    # Iterate over the string
    while i < len(result) - 3:
        # Check each operation
        for pattern, insert_char in operations:
            # If the pattern matches, insert the character
            if result[i:i+4] == list(pattern):
                result.insert(i+4, insert_char)
                i += 4  # Move past the inserted character
                break
        else:
            i += 1  # Move to the next character if no pattern matched
    
    # Join the list back into a string
    return ''.join(result)

# Given string
s = "EABCCCDEAAABCD"
# Apply the T operation
new_string = apply_t_operation(s)
print(new_string)