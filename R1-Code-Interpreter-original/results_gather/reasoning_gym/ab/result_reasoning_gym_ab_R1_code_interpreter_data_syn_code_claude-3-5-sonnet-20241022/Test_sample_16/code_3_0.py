def apply_rules(program):
    tokens = program.split()
    steps = [program]
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Check if we have a pattern where '#' faces '#'
            if tokens[i].endswith('#') and tokens[i+1].startswith('#'):
                pair = tokens[i] + " " + tokens[i+1]
                
                # Apply rules
                if pair == "A# #A" or pair == "B# #B":
                    # Remove both tokens
                    tokens.pop(i)
                    tokens.pop(i)
                    changed = True
                elif pair == "A# #B":
                    tokens[i:i+2] = ["#B", "A#"]
                    changed = True
                elif pair == "B# #A":
                    tokens[i:i+2] = ["#A", "B#"]
                    changed = True
                
                if changed:
                    steps.append(" ".join(tokens))
                    break
            i += 1
        
        if not changed:
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    
    return " ".join(tokens)

initial = "#A B# #A B# #B B# B# A# #B B#"
print(f"Final: {apply_rules(initial)}")