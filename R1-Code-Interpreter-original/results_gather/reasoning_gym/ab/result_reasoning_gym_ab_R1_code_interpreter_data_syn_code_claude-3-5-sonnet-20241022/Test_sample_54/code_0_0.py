def apply_rules(tokens):
    # Convert string to list for easier manipulation
    if isinstance(tokens, str):
        tokens = tokens.split()
    
    # Keep applying rules until no changes can be made
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair of tokens
            pair = tokens[i:i+2]
            
            # Apply rules
            if pair == ['A#', '#A']:
                # Remove both tokens
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                i -= 1  # Move back one position
            elif pair == ['A#', '#B']:
                # Replace with #B A#
                tokens[i:i+2] = ['#B', 'A#']
                changed = True
            elif pair == ['B#', '#A']:
                # Replace with #A B#
                tokens[i:i+2] = ['#A', 'B#']
                changed = True
            elif pair == ['B#', '#B']:
                # Remove both tokens
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                i -= 1  # Move back one position
            else:
                i += 1
                
        if not changed:
            break
    
    return tokens

# Initial program
program = "A# #A A# A# #A #B A# A# #B #B"
result = apply_rules(program)
print(' '.join(result))