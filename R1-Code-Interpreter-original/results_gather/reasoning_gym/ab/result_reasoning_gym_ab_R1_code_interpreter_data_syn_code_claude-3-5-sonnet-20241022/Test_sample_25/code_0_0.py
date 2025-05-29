def apply_rules(program):
    # Convert program string to list of tokens
    tokens = program.split()
    
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
                break
                
            # Rule 2: A# #B → #B A#
            elif pair == "A# #B":
                tokens[i] = "#B"
                tokens[i+1] = "A#"
                changed = True
                break
                
            # Rule 3: B# #A → #A B#
            elif pair == "B# #A":
                tokens[i] = "#A"
                tokens[i+1] = "B#"
                changed = True
                break
                
            # Rule 4: B# #B → nothing
            elif pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
                
            i += 1
            
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "A# #A A# B# #A #B #A #B #B A#"
final_state = apply_rules(program)
print(final_state)