def apply_rules(tokens):
    # Convert string to list for easier manipulation
    if isinstance(tokens, str):
        tokens = tokens.split()
    
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair of tokens
            pair = tokens[i:i+2]
            
            # Rule 1: A# #A -> nothing
            if pair == ['A#', '#A']:
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                continue
                
            # Rule 2: A# #B -> #B A#
            elif pair == ['A#', '#B']:
                tokens[i], tokens[i+1] = tokens[i+1], tokens[i]
                changed = True
                
            # Rule 3: B# #A -> #A B#
            elif pair == ['B#', '#A']:
                tokens[i], tokens[i+1] = tokens[i+1], tokens[i]
                changed = True
                
            # Rule 4: B# #B -> nothing
            elif pair == ['B#', '#B']:
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                continue
                
            i += 1
    
    return tokens

# Initial program
program = "#B #B #B A# #B #A A# #A A# #A"
result = apply_rules(program)
print(' '.join(result))