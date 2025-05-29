def apply_rules(program):
    # Convert program string to list of tokens
    tokens = program.split()
    step = 1
    
    # Keep applying rules until no changes can be made
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i+1]
            
            # Rule 1: A# #A → nothing
            if pair == "A# #A":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                print(f"Step {step}: Applied Rule 1 (A# #A → nothing)")
                print(" ".join(tokens))
                step += 1
                break
                
            # Rule 2: A# #B → #B A#
            elif pair == "A# #B":
                tokens[i] = "#B"
                tokens[i+1] = "A#"
                changed = True
                print(f"Step {step}: Applied Rule 2 (A# #B → #B A#)")
                print(" ".join(tokens))
                step += 1
                break
                
            # Rule 3: B# #A → #A B#
            elif pair == "B# #A":
                tokens[i] = "#A"
                tokens[i+1] = "B#"
                changed = True
                print(f"Step {step}: Applied Rule 3 (B# #A → #A B#)")
                print(" ".join(tokens))
                step += 1
                break
                
            # Rule 4: B# #B → nothing
            elif pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                print(f"Step {step}: Applied Rule 4 (B# #B → nothing)")
                print(" ".join(tokens))
                step += 1
                break
                
            i += 1
            
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "A# #A #B #B B# #A A# B# B# #B"
final_state = apply_rules(program)
print("\nFinal state:", final_state)