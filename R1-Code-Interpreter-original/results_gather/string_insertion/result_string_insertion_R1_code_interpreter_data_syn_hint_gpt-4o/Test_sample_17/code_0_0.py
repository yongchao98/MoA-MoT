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
    while i <= len(s_list) - 4:
        # Check each operation
        for pattern, insert_char in operations:
            # If the pattern is found, insert the character
            if s_list[i:i+4] == list(pattern):
                s_list.insert(i+4, insert_char)
                i += 4  # Move past the inserted character
                break
        i += 1
    
    # Convert the list back to a string
    return ''.join(s_list)

# Given string
initial_string = "ADEABADEABEDEABB"
# Apply the T operation
result_string = apply_t_operation(initial_string)
print(result_string)