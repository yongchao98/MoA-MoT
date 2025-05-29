def apply_rules(program):
    tokens = program.split()
    print("Initial:", " ".join(tokens))
    
    while True:
        i = 0
        changed = False
        while i < len(tokens) - 1:
            pair = f"{tokens[i]} {tokens[i+1]}"
            
            # Apply rules
            if pair in ["A# #A", "B# #B"]:
                tokens = tokens[:i] + tokens[i+2:]
                changed = True
                print("Current:", " ".join(tokens))
                break
            elif pair == "A# #B":
                tokens = tokens[:i] + ["#B", "A#"] + tokens[i+2:]
                changed = True
                print("Current:", " ".join(tokens))
                break
            elif pair == "B# #A":
                tokens = tokens[:i] + ["#A", "B#"] + tokens[i+2:]
                changed = True
                print("Current:", " ".join(tokens))
                break
            i += 1
            
        if not changed:
            break
    
    return " ".join(tokens)

program = "#B #B #A #A #B A# #B A# A# B#"
final = apply_rules(program)
print("\nFinal:", final)