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
            if pair == "A# #A" or pair == "B# #B":
                # Remove both tokens
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            elif pair == "A# #B":
                # Swap positions
                tokens[i], tokens[i+1] = "#B", "A#"
                changed = True
                i += 1
            elif pair == "B# #A":
                # Swap positions
                tokens[i], tokens[i+1] = "#A", "B#"
                changed = True
                i += 1
            else:
                i += 1
        
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "B# B# A# #A #A #A A# #A A# #A"
result = apply_rules(program)
print(result)