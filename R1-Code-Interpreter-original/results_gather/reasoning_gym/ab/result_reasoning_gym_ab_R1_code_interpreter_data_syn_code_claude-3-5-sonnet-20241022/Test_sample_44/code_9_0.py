def apply_rules(program):
    tokens = program.split()
    print("Initial:", ' '.join(tokens))
    
    step = 1
    while True:
        changed = False
        
        # First, handle all adjacent transformations
        i = 0
        while i < len(tokens) - 1:
            # Handle cancellations
            if (tokens[i] == "A#" and tokens[i+1] == "#A") or \
               (tokens[i] == "B#" and tokens[i+1] == "#B"):
                print(f"Step {step}: Cancel {tokens[i]} {tokens[i+1]}")
                tokens.pop(i)
                tokens.pop(i)
                print("After:", ' '.join(tokens))
                changed = True
                break
            
            # Handle B# #A transformations (most important)
            if tokens[i] == "B#" and tokens[i+1] == "#A":
                print(f"Step {step}: Transform B# #A to #A B#")
                tokens[i] = "#A"
                tokens[i+1] = "B#"
                print("After:", ' '.join(tokens))
                changed = True
                break
            
            # Handle A# #B transformations
            if tokens[i] == "A#" and tokens[i+1] == "#B":
                print(f"Step {step}: Transform A# #B to #B A#")
                tokens[i] = "#B"
                tokens[i+1] = "A#"
                print("After:", ' '.join(tokens))
                changed = True
                break
            
            i += 1
        
        # If no adjacent transformations, look for B# that can move left
        if not changed:
            for i in range(len(tokens)-2, -1, -1):
                if tokens[i] == "B#":
                    # Try to move B# left if possible
                    for j in range(i-1, -1, -1):
                        if tokens[j] == "#A":
                            # We can move B# left
                            tokens.pop(i)
                            tokens.insert(j+1, "B#")
                            print(f"Step {step}: Move B# left")
                            print("After:", ' '.join(tokens))
                            changed = True
                            break
                    if changed:
                        break
        
        # If still no change, look for #A that can move right
        if not changed:
            for i in range(len(tokens)-1):
                if tokens[i] == "#A":
                    # Try to move #A right if possible
                    for j in range(i+1, len(tokens)):
                        if tokens[j] == "B#":
                            # We can move #A right
                            tokens.pop(i)
                            tokens.insert(j, "#A")
                            print(f"Step {step}: Move #A right")
                            print("After:", ' '.join(tokens))
                            changed = True
                            break
                    if changed:
                        break
        
        if not changed:
            break
        step += 1
    
    return ' '.join(tokens)

program = "#A #B #B B# #A A# #B B# #A #A"
final_state = apply_rules(program)
print("\nFinal state:", final_state)