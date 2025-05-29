def apply_rules(program):
    # Convert string to list for easier manipulation
    tokens = program.split()
    
    # Keep applying rules until no changes can be made
    while True:
        changed = False
        for i in range(len(tokens) - 1):
            # Check each pair of adjacent tokens
            pair = tokens[i] + " " + tokens[i+1]
            
            # Apply rules
            if pair in ["A# #A", "B# #B"]:
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
program = "#B A# #B A# B# #A B# B# A# #B"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")