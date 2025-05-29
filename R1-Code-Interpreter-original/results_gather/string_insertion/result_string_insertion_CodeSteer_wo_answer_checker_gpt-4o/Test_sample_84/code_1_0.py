def apply_t_operation(s):
    i = 0
    result = ""
    while i < len(s):
        if s[i:i+4] == "ABCD":
            result += "ABCD" + "A"
            i += 4
        elif s[i:i+4] == "BCDE":
            result += "BCDE" + "B"
            i += 4
        elif s[i:i+4] == "CDEA":
            result += "CDEA" + "C"
            i += 4
        elif s[i:i+4] == "DEAB":
            result += "DEAB" + "D"
            i += 4
        elif s[i:i+4] == "EABC":
            result += "EABC" + "E"
            i += 4
        else:
            result += s[i]
            i += 1
    return f"<<<{result}>>>"

# Given string
input_string = "BBCDECBCDEECDEA"
# Apply the T operation
output_string = apply_t_operation(input_string)
# Print the result
print(output_string)