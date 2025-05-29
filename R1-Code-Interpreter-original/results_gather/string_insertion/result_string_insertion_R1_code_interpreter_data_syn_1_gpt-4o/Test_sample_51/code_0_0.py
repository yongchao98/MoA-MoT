def apply_t_operation(s):
    # Define the rules for the T operation
    rules = [
        ("ABCD", "A"),
        ("BCDE", "B"),
        ("CDEA", "C"),
        ("DEAB", "D"),
        ("EABC", "E")
    ]
    
    # Initialize the result string
    result = []
    i = 0
    while i < len(s):
        # Check each rule
        for pattern, insert_char in rules:
            if s[i:i+4] == pattern:
                # If pattern matches, append the pattern and the insert character
                result.append(s[i:i+4] + insert_char)
                i += 4  # Move past the pattern
                break
        else:
            # If no pattern matches, just append the current character
            result.append(s[i])
            i += 1
    
    # Join the result list into a final string
    return ''.join(result)

# Given string
initial_string = "EDDEABDABCDBBCDEC"

# Apply the T operation
new_string = apply_t_operation(initial_string)

# Output the result
print(new_string)