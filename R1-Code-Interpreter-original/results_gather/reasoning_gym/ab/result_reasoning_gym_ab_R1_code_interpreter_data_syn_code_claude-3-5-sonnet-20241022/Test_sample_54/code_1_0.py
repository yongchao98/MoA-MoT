def apply_rules(tokens):
    if isinstance(tokens, str):
        tokens = tokens.split()
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Only apply rules when '#' symbols are adjacent
            current = tokens[i]
            next_token = tokens[i+1]
            
            # Check if '#' symbols are facing each other
            if current.endswith('#') and next_token.startswith('#'):
                if current == 'A#' and next_token == '#A':
                    tokens.pop(i)
                    tokens.pop(i)
                    changed = True
                    i -= 1 if i > 0 else 0
                elif current == 'A#' and next_token == '#B':
                    tokens[i:i+2] = ['#B', 'A#']
                    changed = True
                elif current == 'B#' and next_token == '#A':
                    tokens[i:i+2] = ['#A', 'B#']
                    changed = True
                elif current == 'B#' and next_token == '#B':
                    tokens.pop(i)
                    tokens.pop(i)
                    changed = True
                    i -= 1 if i > 0 else 0
            i += 1
        
        if not changed:
            break
    
    return tokens

program = "A# #A A# A# #A #B A# A# #B #B"
result = apply_rules(program)
print(' '.join(result))