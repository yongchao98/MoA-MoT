def apply_rules(program):
    # Convert string to list for easier manipulation
    tokens = program.split()
    
    # Keep applying rules until no changes can be made
    iteration = 0
    max_iterations = 1000  # Safety limit
    
    while iteration < max_iterations:
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
                i = 0  # Start from beginning after each change
                continue
            elif pair == "A# #B":
                # Swap positions
                tokens[i], tokens[i+1] = "#B", "A#"
                changed = True
                i = 0  # Start from beginning after each change
                continue
            elif pair == "B# #A":
                # Swap positions
                tokens[i], tokens[i+1] = "#A", "B#"
                changed = True
                i = 0  # Start from beginning after each change
                continue
            i += 1
        
        if not changed:
            break
            
        iteration += 1
    
    return " ".join(tokens)

# Initial program
program = "#A #B #B B# #A A# #B B# #A #A"
final_state = apply_rules(program)
print(final_state)