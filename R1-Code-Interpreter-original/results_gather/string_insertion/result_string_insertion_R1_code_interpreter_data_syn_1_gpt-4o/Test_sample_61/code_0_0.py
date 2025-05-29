def apply_t_operation(s):
    insertions = []
    i = 0
    while i <= len(s) - 4:
        if s[i:i+4] == "ABCD":
            insertions.append((i+4, 'A'))
        elif s[i:i+4] == "BCDE":
            insertions.append((i+4, 'B'))
        elif s[i:i+4] == "CDEA":
            insertions.append((i+4, 'C'))
        elif s[i:i+4] == "DEAB":
            insertions.append((i+4, 'D'))
        elif s[i:i+4] == "EABC":
            insertions.append((i+4, 'E'))
        i += 1

    # Insert characters from left to right
    for index, char in sorted(insertions):
        s = s[:index] + char + s[index:]
    
    return s

# Given string
s = "DEABDCDEADABCDA"
# Apply the T operation
result = apply_t_operation(s)
print(result)