def apply_rules(program):
    tokens = program.split()
    
    while True:
        original = tokens.copy()
        i = 0
        while i < len(tokens) - 1:
            pair = f"{tokens[i]} {tokens[i+1]}"
            
            if pair == "A# #A" or pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                continue
            elif pair == "A# #B":
                tokens[i] = "#B"
                tokens[i+1] = "A#"
            elif pair == "B# #A":
                tokens[i] = "#A"
                tokens[i+1] = "B#"
            i += 1
        
        if tokens == original:
            break
    
    return " ".join(tokens)

program = "#A A# #A B# B# #A A# B# A# B#"
result = apply_rules(program)
print(result)