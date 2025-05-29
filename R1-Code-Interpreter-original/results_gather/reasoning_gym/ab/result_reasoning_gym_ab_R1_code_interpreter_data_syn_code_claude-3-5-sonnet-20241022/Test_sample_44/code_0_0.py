def apply_rules(program):
    # Convert string to list for easier manipulation
    tokens = program.split()
    
    # Keep applying rules until no changes can be made
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair
            pair = tokens[i] + " " + tokens[i+1]
            
            # Apply rules
            if pair == "A# #A" or pair == "B# #B":
                # Remove both tokens
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            elif pair == "A# #B":
                # Swap positions
                tokens[i], tokens[i+1] = "#B", "A#"
                changed = True
                break
            elif pair == "B# #A":
                # Swap positions
                tokens[i], tokens[i+1] = "#A", "B#"
                changed = True
                break
            i += 1
        
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "#A #B #B B# #A A# #B B# #A #A"
final_state = apply_rules(program)
print(final_state)