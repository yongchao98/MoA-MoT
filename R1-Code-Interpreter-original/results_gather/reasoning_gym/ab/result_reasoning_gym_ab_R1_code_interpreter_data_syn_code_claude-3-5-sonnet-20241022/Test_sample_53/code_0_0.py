def apply_rules(program):
    # Convert string to list for easier manipulation
    tokens = program.split()
    
    # Keep applying rules until no changes can be made
    while True:
        changed = False
        for i in range(len(tokens) - 1):
            # Skip if we've removed this token
            if tokens[i] is None or tokens[i+1] is None:
                continue
                
            # Check each rule
            pair = f"{tokens[i]} {tokens[i+1]}"
            
            if pair == "A# #A" or pair == "B# #B":
                # Remove both tokens
                tokens[i] = None
                tokens[i+1] = None
                changed = True
                
            elif pair == "A# #B":
                # Replace with "#B A#"
                tokens[i] = "#B"
                tokens[i+1] = "A#"
                changed = True
                
            elif pair == "B# #A":
                # Replace with "#A B#"
                tokens[i] = "#A"
                tokens[i+1] = "B#"
                changed = True
        
        # Remove None values
        tokens = [t for t in tokens if t is not None]
        
        # If no changes were made, we're done
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "B# A# #A B# B# #B #A B# #B B#"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")