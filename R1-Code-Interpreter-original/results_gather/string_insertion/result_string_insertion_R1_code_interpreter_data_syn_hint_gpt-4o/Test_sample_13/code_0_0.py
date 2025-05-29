def apply_t_operation(s):
    # Define the substrings and their corresponding insertions
    operations = [
        ("ABCD", "A"),
        ("BCDE", "B"),
        ("CDEA", "C"),
        ("DEAB", "D"),
        ("EABC", "E")
    ]
    
    # Initialize the result as a list of characters for easier manipulation
    result = list(s)
    i = 0
    
    # Iterate through the string
    while i < len(result) - 3:
        # Check each operation
        for substring, char_to_insert in operations:
            # If the substring is found, insert the character
            if result[i:i+4] == list(substring):
                result.insert(i+4, char_to_insert)
                i += 4  # Move past the inserted character
                break
        else:
            i += 1  # Move to the next character if no insertion was made
    
    # Join the list back into a string
    return ''.join(result)

# Given string
s = "DEABDDEABADEABBA"
# Apply the T operation
resulting_string = apply_t_operation(s)
print(resulting_string)