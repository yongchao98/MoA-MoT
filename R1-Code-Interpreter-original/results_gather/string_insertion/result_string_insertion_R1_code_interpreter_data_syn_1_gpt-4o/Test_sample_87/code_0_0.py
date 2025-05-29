def apply_t_operation(s):
    # Define the substrings and their corresponding insertions
    operations = [
        ("ABCD", "A"),
        ("BCDE", "B"),
        ("CDEA", "C"),
        ("DEAB", "D"),
        ("EABC", "E")
    ]
    
    # List to keep track of insertions
    insertions = []
    
    # Iterate over the string to find all substrings
    for i in range(len(s) - 3):
        for sub, char in operations:
            if s[i:i+4] == sub:
                insertions.append((i+4, char))
    
    # Sort insertions by index to ensure left-to-right order
    insertions.sort()
    
    # Apply insertions
    result = []
    last_index = 0
    for index, char in insertions:
        result.append(s[last_index:index])
        result.append(char)
        last_index = index
    
    # Append the remaining part of the string
    result.append(s[last_index:])
    
    return ''.join(result)

# Given string
s = "ADABCDCCDEADABCD"
# Apply the T operation
result = apply_t_operation(s)
print(result)