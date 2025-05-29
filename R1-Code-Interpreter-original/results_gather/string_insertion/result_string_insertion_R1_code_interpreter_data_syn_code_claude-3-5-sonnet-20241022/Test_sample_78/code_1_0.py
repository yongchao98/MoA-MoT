def apply_T_with_steps(s):
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    result = list(s)
    i = 0
    
    while i <= len(result) - 4:
        substring = ''.join(result[i:i+4])
        if substring in patterns:
            print(f"Found {substring} at position {i}")
            print(f"Before: {''.join(result)}")
            result.insert(i + 4, patterns[substring])
            print(f"After:  {''.join(result)}\n")
            i += 5
        else:
            i += 1
    
    return ''.join(result)

input_string = "ABCDBABCDCCDEA"
result = apply_T_with_steps(input_string)
print("Final result:", result)