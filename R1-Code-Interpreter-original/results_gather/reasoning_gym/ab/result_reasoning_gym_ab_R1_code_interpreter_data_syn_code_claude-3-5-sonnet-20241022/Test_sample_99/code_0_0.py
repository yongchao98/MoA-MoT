def can_reduce(tokens, i):
    if i >= len(tokens) - 1:
        return False, None
    
    current = tokens[i]
    next_token = tokens[i + 1]
    
    # Check for disappearing rules
    if (current == "A#" and next_token == "#A") or \
       (current == "B#" and next_token == "#B"):
        return True, ("delete", i)
    
    # Check for swapping rules
    if (current == "A#" and next_token == "#B") or \
       (current == "B#" and next_token == "#A"):
        return True, ("swap", i)
    
    return False, None

def process_program(program):
    # Split into tokens
    tokens = program.split()
    
    while True:
        made_change = False
        i = 0
        while i < len(tokens):
            can_change, action = can_reduce(tokens, i)
            if can_change:
                made_change = True
                if action[0] == "delete":
                    # Remove both tokens
                    tokens.pop(i)
                    tokens.pop(i)
                else:  # swap
                    tokens[i], tokens[i + 1] = tokens[i + 1], tokens[i]
                    i += 1
            else:
                i += 1
        
        if not made_change:
            break
    
    return " ".join(tokens)

# Initial program
program = "#A #A A# #B #A A# #B A# #A A#"
result = process_program(program)
print(result)