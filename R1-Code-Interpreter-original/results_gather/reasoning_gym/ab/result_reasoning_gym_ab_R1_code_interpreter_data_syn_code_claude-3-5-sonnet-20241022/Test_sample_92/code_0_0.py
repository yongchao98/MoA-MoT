def apply_rules(program):
    # Convert string to list for easier manipulation
    tokens = program.split()
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair of tokens
            pair = tokens[i] + " " + tokens[i+1]
            
            # Apply rules
            if pair in ["A# #A", "B# #B"]:
                # Remove both tokens
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                i -= 1  # Adjust index after removal
            elif pair == "A# #B":
                # Replace with "#B A#"
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
            elif pair == "B# #A":
                # Replace with "#A B#"
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
            i += 1
            
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "#A A# B# A# #B #A #A A# #A B#"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")