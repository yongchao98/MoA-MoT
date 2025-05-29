def apply_rules(tokens):
    tokens = tokens.split()
    step = 1
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Check if we have a B# #A pattern (ignoring spaces)
            if tokens[i].endswith('#') and tokens[i+1].startswith('#'):
                pair = f"{tokens[i]} {tokens[i+1]}"
                
                if pair == "A# #A" or pair == "B# #B":
                    tokens.pop(i)
                    tokens.pop(i)
                    changed = True
                    break
                elif pair == "A# #B":
                    tokens[i:i+2] = ["#B", "A#"]
                    changed = True
                    break
                elif pair == "B# #A":
                    tokens[i:i+2] = ["#A", "B#"]
                    changed = True
                    break
            i += 1
            
        if changed:
            print(f"Step {step}: {' '.join(tokens)}")
            step += 1
        else:
            break
    
    return " ".join(tokens)

program = "#B #A #B A# B# A# A# A#"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")