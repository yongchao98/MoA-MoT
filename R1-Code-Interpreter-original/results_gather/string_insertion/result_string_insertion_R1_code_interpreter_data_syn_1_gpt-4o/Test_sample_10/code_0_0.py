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
        for pattern, insert_char in operations:
            # Check if the current substring matches the pattern
            if result[i:i+4] == list(pattern):
                # Insert the character after the pattern
                result.insert(i+4, insert_char)
                # Move the index past the inserted character
                i += 4
                break
        i += 1
    
    # Join the list back into a string
    return ''.join(result)

# Given string
s = "BABCDDBCDECEABC"
# Apply the T operation
result = apply_t_operation(s)
print(result)