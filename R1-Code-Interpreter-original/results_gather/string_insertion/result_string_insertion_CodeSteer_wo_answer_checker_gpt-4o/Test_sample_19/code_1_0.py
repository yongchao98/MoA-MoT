def apply_t_operation(s):
    i = 0
    while i < len(s) - 3:
        if s[i:i+4] == "ABCD":
            s = s[:i+4] + "A" + s[i+4:]
            i += 5  # Move past the inserted character
        elif s[i:i+4] == "BCDE":
            s = s[:i+4] + "B" + s[i+4:]
            i += 5
        elif s[i:i+4] == "CDEA":
            s = s[:i+4] + "C" + s[i+4:]
            i += 5
        elif s[i:i+4] == "DEAB":
            s = s[:i+4] + "D" + s[i+4:]
            i += 5
        elif s[i:i+4] == "EABC":
            s = s[:i+4] + "E" + s[i+4:]
            i += 5
        else:
            i += 1
    return f"<<<{s}>>>"

# Given string
input_string = "BCDEEEABCCEABCC"
# Apply the T operation
result = apply_t_operation(input_string)
# Print the result
print(result)