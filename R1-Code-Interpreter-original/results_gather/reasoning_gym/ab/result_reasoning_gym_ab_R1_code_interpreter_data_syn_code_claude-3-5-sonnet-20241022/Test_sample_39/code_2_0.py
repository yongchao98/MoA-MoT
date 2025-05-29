def apply_rules(program):
    tokens = program.split()
    
    def find_match(tokens):
        for i in range(len(tokens) - 1):
            pair = f"{tokens[i]} {tokens[i+1]}"
            if pair in ["A# #A", "B# #B"]:
                return i, "remove"
            elif pair == "A# #B":
                return i, "swap1"
            elif pair == "B# #A":
                return i, "swap2"
        return -1, None

    while True:
        pos, action = find_match(tokens)
        if pos == -1:
            break
            
        if action == "remove":
            tokens.pop(pos)
            tokens.pop(pos)
        elif action == "swap1":
            tokens[pos] = "#B"
            tokens[pos+1] = "A#"
        elif action == "swap2":
            tokens[pos] = "#A"
            tokens[pos+1] = "B#"
    
    return " ".join(tokens)

program = "#A A# #A B# B# #A A# B# A# B#"
result = apply_rules(program)
print(result)