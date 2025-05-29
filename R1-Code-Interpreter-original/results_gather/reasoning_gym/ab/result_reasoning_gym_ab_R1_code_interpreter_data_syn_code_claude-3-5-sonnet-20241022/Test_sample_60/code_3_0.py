def apply_rules(tokens):
    tokens = tokens.split()
    step = 1
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair of tokens
            t1, t2 = tokens[i], tokens[i+1]
            
            # Check if this is a pattern where '#' symbols face each other
            if (t1.endswith('#') and t2.startswith('#')) or \
               (t1 == "B#" and t2.startswith('#')) or \
               (t1.endswith('#') and t2 == "#A"):
                
                if t1 == "A#" and t2 == "#A":
                    tokens.pop(i)
                    tokens.pop(i)
                    changed = True
                    break
                elif t1 == "B#" and t2 == "#B":
                    tokens.pop(i)
                    tokens.pop(i)
                    changed = True
                    break
                elif t1 == "A#" and t2 == "#B":
                    tokens[i:i+2] = ["#B", "A#"]
                    changed = True
                    break
                elif t1 == "B#" and t2 == "#A":
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