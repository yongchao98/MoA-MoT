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
            i += 1
        
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "#B #A #B #B A# B# B# #B A# #A"
result = apply_rules(program)
print(result)