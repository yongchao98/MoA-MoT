def apply_rules(tokens):
    # Convert string to list for easier manipulation
    tokens = tokens.split()
    
    # Keep applying rules until no changes can be made
    while True:
        changed = False
        for i in range(len(tokens) - 1):
            # Check each pair of adjacent tokens
            pair = tokens[i] + " " + tokens[i+1]
            
            # Apply rules
            if pair == "A# #A" or pair == "B# #B":
                # Remove both tokens
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            elif pair == "A# #B":
                # Replace with "#B A#"
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
                break
            elif pair == "B# #A":
                # Replace with "#A B#"
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
                break
        
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "#A B# A# A# #A #A #A #A #B A#"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")