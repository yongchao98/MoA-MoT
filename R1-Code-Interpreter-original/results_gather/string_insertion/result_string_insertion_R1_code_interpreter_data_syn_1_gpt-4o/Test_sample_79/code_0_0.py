def apply_t_operation(s):
    i = 0
    result = []
    while i < len(s):
        if s[i:i+4] == "ABCD":
            result.append(s[i:i+4] + "A")
            i += 4
        else:
            result.append(s[i])
            i += 1
    return ''.join(result)

initial_string = "DABCDBABCDBABCDB"
final_string = apply_t_operation(initial_string)
print(final_string)