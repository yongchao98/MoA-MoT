def apply_rules(tokens):
    # Returns (new_tokens, changed)
    # where changed indicates if any rule was applied
    for i in range(len(tokens) - 1):
        pair = tokens[i:i+2]
        
        # Rule 1: A# #A -> nothing
        if pair == ['A#', '#A']:
            return tokens[:i] + tokens[i+2:], True
        
        # Rule 2: A# #B -> #B A#
        elif pair == ['A#', '#B']:
            return tokens[:i] + ['#B', 'A#'] + tokens[i+2:], True
            
        # Rule 3: B# #A -> #A B#
        elif pair == ['B#', '#A']:
            return tokens[:i] + ['#A', 'B#'] + tokens[i+2:], True
            
        # Rule 4: B# #B -> nothing
        elif pair == ['B#', '#B']:
            return tokens[:i] + tokens[i+2:], True
    
    return tokens, False

# Initial program
program = ['A#', '#A', '#B', 'B#', '#B', '#A', '#A', 'A#', '#B', 'B#']

# Keep applying rules until no changes occur
while True:
    new_program, changed = apply_rules(program)
    if not changed:
        break
    program = new_program

# Print the final state
print(' '.join(program))