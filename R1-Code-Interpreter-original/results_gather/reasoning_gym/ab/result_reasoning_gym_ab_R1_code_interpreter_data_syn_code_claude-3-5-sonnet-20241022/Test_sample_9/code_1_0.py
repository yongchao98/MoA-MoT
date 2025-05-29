def apply_rules(program):
    tokens = program.split()
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            current = tokens[i]
            next_token = tokens[i + 1]
            
            # Special case for #A A# (equivalent to A# #A)
            if (current == "#A" and next_token == "A#"):
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            
            # Original rules
            elif current == "A#" and next_token == "#A":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            elif current == "A#" and next_token == "#B":
                tokens[i] = "#B"
                tokens[i + 1] = "A#"
                changed = True
                break
            elif current == "B#" and next_token == "#A":
                tokens[i] = "#A"
                tokens[i + 1] = "B#"
                changed = True
                break
            elif current == "B#" and next_token == "#B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            
            i += 1
        
        if not changed:
            break
    
    return " ".join(tokens)

program = "#B #A #A #A A# B#"
final_state = apply_rules(program)
print(final_state)