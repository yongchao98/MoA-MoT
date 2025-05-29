def apply_t_operation(s):
    insertions = []
    i = 0
    while i < len(s) - 3:
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

    # Apply insertions from left to right
    result = []
    last_index = 0
    for index, char in insertions:
        result.append(s[last_index:index])
        result.append(char)
        last_index = index
    result.append(s[last_index:])

    return ''.join(result)

# Given string
s = "DACDEAACDEADDEAB"
# Apply the T operation
result = apply_t_operation(s)
print(result)