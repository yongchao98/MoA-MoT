def process_program(program):
    # Convert string to list of tokens
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
                break
            elif pair == "A# #B":
                # Swap tokens
                tokens[i], tokens[i+1] = "#B", "A#"
                changed = True
                break
            elif pair == "B# #A":
                # Swap tokens
                tokens[i], tokens[i+1] = "#A", "B#"
                changed = True
                break
            i += 1
        
        # If no changes were made, we're done
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "B# A# #A #A #B A# #B #A #A #B"
print(process_program(program))