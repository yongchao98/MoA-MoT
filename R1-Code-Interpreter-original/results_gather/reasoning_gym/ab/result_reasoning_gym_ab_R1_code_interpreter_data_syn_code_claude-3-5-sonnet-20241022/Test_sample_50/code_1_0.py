def apply_rules(program):
    tokens = program.split()
    steps = [" ".join(tokens)]
    
    while True:
        changed = False
        for i in range(len(tokens) - 1):
            pair = f"{tokens[i]} {tokens[i+1]}"
            
            # Rule 1: A# #A -> nothing
            if pair == "A# #A":
                tokens = tokens[:i] + tokens[i+2:]
                changed = True
                break
                
            # Rule 2: A# #B -> #B A#
            elif pair == "A# #B":
                tokens = tokens[:i] + ["#B", "A#"] + tokens[i+2:]
                changed = True
                break
                
            # Rule 3: B# #A -> #A B#
            elif pair == "B# #A":
                tokens = tokens[:i] + ["#A", "B#"] + tokens[i+2:]
                changed = True
                break
                
            # Rule 4: B# #B -> nothing
            elif pair == "B# #B":
                tokens = tokens[:i] + tokens[i+2:]
                changed = True
                break
        
        if changed:
            steps.append(" ".join(tokens))
        else:
            break
    
    print("\nTransformation steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return " ".join(tokens)

program = "B# #B #A B# B# #A B# A# #B #B"
final_state = apply_rules(program)
print("\nFinal state:", final_state)