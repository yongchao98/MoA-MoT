def apply_rules(tokens):
    # Convert to list for easier manipulation
    tokens = tokens.split()
    
    # Keep applying rules until no changes can be made
    while True:
        changed = False
        for i in range(len(tokens) - 1):
            # Skip if we've deleted this token
            if tokens[i] is None or tokens[i+1] is None:
                continue
                
            # Check for facing '#' symbols
            if tokens[i].endswith('#') and tokens[i+1].startswith('#'):
                # Apply rules
                if tokens[i] == 'A#' and tokens[i+1] == '#A':
                    tokens[i] = tokens[i+1] = None
                    changed = True
                elif tokens[i] == 'A#' and tokens[i+1] == '#B':
                    tokens[i], tokens[i+1] = '#B', 'A#'
                    changed = True
                elif tokens[i] == 'B#' and tokens[i+1] == '#A':
                    tokens[i], tokens[i+1] = '#A', 'B#'
                    changed = True
                elif tokens[i] == 'B#' and tokens[i+1] == '#B':
                    tokens[i] = tokens[i+1] = None
                    changed = True
        
        # Remove None values
        tokens = [t for t in tokens if t is not None]
        
        # Print current state
        print("Current state:", ' '.join(tokens))
        
        # If no changes were made, we're done
        if not changed:
            break
    
    return ' '.join(tokens)

# Initial program
program = "B# #A #B A# B# B# #A #B B# #A"
result = apply_rules(program)
print("\nFinal result:", result)