def apply_t_operation(s):
    i = 0
    while i < len(s) - 3:
        if s[i:i+4] == "CDEA":
            s = s[:i+4] + "C" + s[i+4:]
            i += 5
        elif s[i:i+4] == "DEAB":
            s = s[:i+4] + "D" + s[i+4:]
            i += 5
        else:
            i += 1
    return s

result = apply_t_operation("ECDEAACDEAEDEAB")
print(result)