def apply_rules(tokens):
    # Returns (new_tokens, changed) where changed indicates if any rule was applied
    i = 0
    while i < len(tokens) - 1:
        # Get current pair of tokens
        pair = tokens[i:i+2]
        
        # Check if '#' symbols face each other
        if ('#' in pair[0] and '#' in pair[1] and 
            (pair[0].index('#') == len(pair[0])-1) != (pair[1].index('#') == len(pair[1])-1)):
            
            # Apply rules
            if pair == ['A#', '#A'] or pair == ['B#', '#B']:
                # Remove both tokens
                tokens = tokens[:i] + tokens[i+2:]
                return tokens, True
            elif pair == ['A#', '#B']:
                # Rewrite as #B A#
                tokens = tokens[:i] + ['#B', 'A#'] + tokens[i+2:]
                return tokens, True
            elif pair == ['B#', '#A']:
                # Rewrite as #A B#
                tokens = tokens[:i] + ['#A', 'B#'] + tokens[i+2:]
                return tokens, True
        i += 1
    return tokens, False

# Initial program
program = ['#B', '#A', '#A', 'B#', '#B', '#B', '#B', '#A', '#B', 'B#']

# Keep applying rules until no more changes
while True:
    program, changed = apply_rules(program)
    if not changed:
        break

print(' '.join(program))