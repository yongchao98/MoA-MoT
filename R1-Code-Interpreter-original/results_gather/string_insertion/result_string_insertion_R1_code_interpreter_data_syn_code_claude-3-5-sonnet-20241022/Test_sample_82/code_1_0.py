def apply_T(s):
    # Define the rules as tuples (pattern, insert_char)
    rules = [
        ('ABCD', 'A'),
        ('BCDE', 'B'),
        ('CDEA', 'C'),
        ('DEAB', 'D'),
        ('EABC', 'E')
    ]
    
    result = list(s)
    insertions = 0
    
    # Print intermediate steps
    for i in range(len(s)):
        substr = s[i:i+4]
        for pattern, insert_char in rules:
            if substr == pattern:
                insert_pos = i + 4 + insertions
                result.insert(insert_pos, insert_char)
                insertions += 1
                print(f"Found {pattern} at position {i}, inserting {insert_char}")
                print(f"Current string: {''.join(result)}")
    
    return ''.join(result)

input_str = "ABCDBDEABAABCD"
result = apply_T(input_str)
print("\nFinal result:", result)