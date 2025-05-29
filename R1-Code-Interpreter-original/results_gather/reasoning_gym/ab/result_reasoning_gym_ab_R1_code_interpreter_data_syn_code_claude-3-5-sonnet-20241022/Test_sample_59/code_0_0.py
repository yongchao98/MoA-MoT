def apply_rules(tokens):
    # Returns (new_tokens, changed) where changed indicates if any rule was applied
    i = 0
    while i < len(tokens) - 1:
        # Get current pair of tokens
        pair = tokens[i:i+2]
        
        # Check if '#' symbols face each other
        if ('#' in pair[0] and '#' in pair[1] and 
            pair[0].index('#') == len(pair[0])-1 and 
            pair[1].index('#') == 0):
            
            # Apply rules
            if pair == ['A#', '#A'] or pair == ['B#', '#B']:
                # Remove both tokens
                tokens = tokens[:i] + tokens[i+2:]
                return tokens, True
            elif pair == ['A#', '#B']:
                # Replace with #B A#
                tokens = tokens[:i] + ['#B', 'A#'] + tokens[i+2:]
                return tokens, True
            elif pair == ['B#', '#A']:
                # Replace with #A B#
                tokens = tokens[:i] + ['#A', 'B#'] + tokens[i+2:]
                return tokens, True
        i += 1
    return tokens, False

# Initial program
program = ['A#', 'B#', 'B#', '#B', '#B', '#B', '#B', '#B', '#A', 'B#']

# Keep applying rules until no more changes
changed = True
while changed:
    program, changed = apply_rules(program)
    print("Current state:", ' '.join(program))

print("\nFinal state:", ' '.join(program))