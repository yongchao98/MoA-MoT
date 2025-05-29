def process_program(program):
    tokens = program.split()
    changed = True
    
    while changed:
        changed = False
        i = 0
        new_tokens = []
        
        while i < len(tokens):
            if i < len(tokens) - 1:
                t1, t2 = tokens[i], tokens[i+1]
                
                if t1 == 'A#' and t2 == '#A':
                    i += 2
                    changed = True
                    continue
                elif t1 == 'B#' and t2 == '#B':
                    i += 2
                    changed = True
                    continue
                elif t1 == 'A#' and t2 == '#B':
                    new_tokens.extend(['#B', 'A#'])
                    i += 2
                    changed = True
                elif t1 == 'B#' and t2 == '#A':
                    new_tokens.extend(['#A', 'B#'])
                    i += 2
                    changed = True
                else:
                    new_tokens.append(t1)
                    i += 1
            else:
                new_tokens.append(tokens[i])
                i += 1
        
        if not changed:
            return ' '.join(tokens)
        tokens = new_tokens
    
    return ' '.join(tokens)

program = "A# B# #A #A #A B# #A B# #B B#"
print(process_program(program))