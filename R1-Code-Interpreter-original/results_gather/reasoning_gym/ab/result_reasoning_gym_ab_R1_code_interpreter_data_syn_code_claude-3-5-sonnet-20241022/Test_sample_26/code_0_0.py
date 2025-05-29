def apply_rules(tokens):
    # Returns (new_tokens, changed)
    # where changed indicates if any rule was applied
    for i in range(len(tokens) - 1):
        # Check if we have a pair of tokens where '#' symbols face each other
        token1 = tokens[i]
        token2 = tokens[i + 1]
        
        if token1.endswith('#') and token2.startswith('#'):
            # Create new list without the matched pair
            new_tokens = tokens[:i]
            
            # Apply the appropriate rule
            if (token1, token2) == ('A#', '#A') or (token1, token2) == ('B#', '#B'):
                # These pairs disappear
                new_tokens.extend(tokens[i + 2:])
            elif (token1, token2) == ('A#', '#B'):
                # A# #B becomes #B A#
                new_tokens.extend(['#B', 'A#'])
                new_tokens.extend(tokens[i + 2:])
            elif (token1, token2) == ('B#', '#A'):
                # B# #A becomes #A B#
                new_tokens.extend(['#A', 'B#'])
                new_tokens.extend(tokens[i + 2:])
            
            return new_tokens, True
    
    return tokens, False

# Initial program
program = ['B#', 'B#', 'B#', 'A#', 'B#', '#B', '#B', '#A', '#A', '#A']

# Keep applying rules until no more changes
while True:
    program, changed = apply_rules(program)
    if not changed:
        break

print(' '.join(program))