def process_program(program):
    # Convert string to list of tokens for easier manipulation
    tokens = program.split()
    
    # Keep processing until no changes occur
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
                i -= 1  # Adjust index after removal
            elif pair == "A# #B":
                # Swap tokens
                tokens[i], tokens[i+1] = "#B", "A#"
                changed = True
            elif pair == "B# #A":
                # Swap tokens
                tokens[i], tokens[i+1] = "#A", "B#"
                changed = True
            i += 1
            
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "B# #A A# #A #A A# #A #A A# #B"
result = process_program(program)
print(result)