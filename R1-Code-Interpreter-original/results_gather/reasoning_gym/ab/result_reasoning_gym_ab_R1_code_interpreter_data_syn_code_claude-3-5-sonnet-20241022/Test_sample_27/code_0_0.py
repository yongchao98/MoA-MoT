def can_apply_rule(tokens, i):
    if i + 1 >= len(tokens):
        return False
    
    # Get the pair of tokens
    t1, t2 = tokens[i], tokens[i+1]
    
    # Check if '#' symbols face each other
    # For t1, '#' should be at the end
    # For t2, '#' should be at the start
    if '#' != t1[-1] or '#' != t2[0]:
        return False
        
    return True

def apply_rule(tokens, i):
    t1, t2 = tokens[i], tokens[i+1]
    
    if t1 == 'A#' and t2 == '#A':
        return [], 2  # Remove both tokens
    elif t1 == 'A#' and t2 == '#B':
        return ['#B', 'A#'], 0  # Swap positions
    elif t1 == 'B#' and t2 == '#A':
        return ['#A', 'B#'], 0  # Swap positions
    elif t1 == 'B#' and t2 == '#B':
        return [], 2  # Remove both tokens
    
    return [t1, t2], 0  # No rule applies

def compute_program(program):
    tokens = program.split()
    
    while True:
        changes_made = False
        i = 0
        
        while i < len(tokens):
            if can_apply_rule(tokens, i):
                new_tokens, skip = apply_rule(tokens, i)
                if skip > 0:
                    tokens[i:i+2] = new_tokens
                    changes_made = True
                    i += 1
                elif new_tokens != tokens[i:i+2]:
                    tokens[i:i+2] = new_tokens
                    changes_made = True
                else:
                    i += 1
            else:
                i += 1
        
        if not changes_made:
            break
    
    return ' '.join(tokens)

# Test the program
program = "B# A# #B A# B# #A A# A# #B #B"
result = compute_program(program)
print(result)