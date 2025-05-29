def apply_t_operation(s):
    # Define the substrings and their corresponding insertions
    operations = {
        "ABCD": "A",
        "BCDE": "B",
        "CDEA": "C",
        "DEAB": "D",
        "EABC": "E"
    }
    
    # Initialize the result as a list of characters for easier manipulation
    result = list(s)
    i = 0
    
    # Process the string from left to right
    while i <= len(result) - 4:
        # Check each substring and perform the corresponding insertion
        for sub, char in operations.items():
            if result[i:i+4] == list(sub):
                result.insert(i+4, char)
                i += 4  # Move past the inserted character
                break
        else:
            i += 1  # Move to the next character if no insertion was made
    
    # Join the list back into a string
    return ''.join(result)

# Given string
s = "EBCDEEBCDEEDEABA"
# Apply the T operation
result = apply_t_operation(s)
print(result)