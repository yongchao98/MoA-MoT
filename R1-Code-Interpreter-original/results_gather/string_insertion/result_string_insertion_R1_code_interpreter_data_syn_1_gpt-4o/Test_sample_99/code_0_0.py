def apply_t_operation(s):
    i = 0
    result = []
    while i < len(s):
        if s[i:i+4] == "ABCD":
            result.extend(s[i:i+4])
            result.append('A')
            i += 4
        elif s[i:i+4] == "BCDE":
            result.extend(s[i:i+4])
            result.append('B')
            i += 4
        elif s[i:i+4] == "CDEA":
            result.extend(s[i:i+4])
            result.append('C')
            i += 4
        elif s[i:i+4] == "DEAB":
            result.extend(s[i:i+4])
            result.append('D')
            i += 4
        elif s[i:i+4] == "EABC":
            result.extend(s[i:i+4])
            result.append('E')
            i += 4
        else:
            result.append(s[i])
            i += 1
    return ''.join(result)

# Given string
s = "CADEABACDEADCDEAED"
# Apply the T operation
new_string = apply_t_operation(s)
print(new_string)