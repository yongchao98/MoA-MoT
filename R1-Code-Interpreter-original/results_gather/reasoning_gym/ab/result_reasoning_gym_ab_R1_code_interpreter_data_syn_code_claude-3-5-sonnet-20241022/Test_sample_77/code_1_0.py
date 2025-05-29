def can_reduce(tokens, pos):
    if pos + 1 >= len(tokens):
        return False
    
    current = tokens[pos]
    next_token = tokens[pos + 1]
    
    # Check if '#' are facing each other
    if current.endswith('#') and next_token.startswith('#'):
        return True
    return False

def apply_rule(tokens, pos):
    current = tokens[pos]
    next_token = tokens[pos + 1]
    
    # Rule: A# #A or B# #B becomes nothing
    if (current == 'A#' and next_token == '#A') or \
       (current == 'B#' and next_token == '#B'):
        return tokens[:pos] + tokens[pos+2:]
    
    # Rule: A# #B becomes #B A#
    elif current == 'A#' and next_token == '#B':
        return tokens[:pos] + ['#B', 'A#'] + tokens[pos+2:]
    
    # Rule: B# #A becomes #A B#
    elif current == 'B#' and next_token == '#A':
        return tokens[:pos] + ['#A', 'B#'] + tokens[pos+2:]
    
    return tokens

def reduce_program(program):
    tokens = program.split()
    step = 1
    
    while True:
        found_reduction = False
        for i in range(len(tokens)):
            if can_reduce(tokens, i):
                new_tokens = apply_rule(tokens, i)
                if new_tokens != tokens:
                    tokens = new_tokens
                    print(f"Step {step}: {' '.join(tokens)}")
                    step += 1
                    found_reduction = True
                    break
        
        if not found_reduction:
            break
    
    return ' '.join(tokens)

# Initial program
program = "#A B# #B #B A# A# B# #B #A A#"
print(f"Initial: {program}")
result = reduce_program(program)
print(f"Final: {result}")