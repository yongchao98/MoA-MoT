def apply_rules(tokens):
    # Convert string to list for easier manipulation
    if isinstance(tokens, str):
        tokens = tokens.split()
    
    # Keep track if any changes were made
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair of tokens
            pair = tokens[i:i+2]
            
            # Check if '#' symbols face each other
            if pair[0].endswith('#') and pair[1].startswith('#'):
                # Apply rules
                if pair == ['A#', '#A'] or pair == ['B#', '#B']:
                    # Remove both tokens
                    tokens.pop(i)
                    tokens.pop(i)
                    changed = True
                elif pair == ['A#', '#B']:
                    # Replace with '#B A#'
                    tokens[i:i+2] = ['#B', 'A#']
                    changed = True
                elif pair == ['B#', '#A']:
                    # Replace with '#A B#'
                    tokens[i:i+2] = ['#A', 'B#']
                    changed = True
                if changed:
                    break
            i += 1
    
    return tokens

# Initial program
program = "A# A# #A #B A# #B #A #B #B #B"
result = apply_rules(program)
print(" ".join(result))