def can_process(tokens, i):
    if i >= len(tokens) - 1:
        return False
    t1, t2 = tokens[i], tokens[i+1]
    return (t1 in ['A#', 'B#']) and (t2 in ['#A', '#B'])

def process_pair(t1, t2):
    if t1 == 'A#' and t2 == '#A':
        return []  # both disappear
    if t1 == 'B#' and t2 == '#B':
        return []  # both disappear
    if t1 == 'A#' and t2 == '#B':
        return ['#B', 'A#']  # swap
    if t1 == 'B#' and t2 == '#A':
        return ['#A', 'B#']  # swap
    return [t1, t2]  # no change

def process_program(program):
    tokens = program.split()
    changed = True
    
    while changed:
        changed = False
        i = 0
        new_tokens = []
        
        while i < len(tokens):
            if i < len(tokens) - 1 and can_process(tokens, i):
                result = process_pair(tokens[i], tokens[i+1])
                new_tokens.extend(result)
                if result != [tokens[i], tokens[i+1]]:
                    changed = True
                i += 2
            else:
                new_tokens.append(tokens[i])
                i += 1
        
        tokens = new_tokens
    
    return ' '.join(tokens)

program = "A# B# #A #A #A B# #A B# #B B#"
print(process_program(program))