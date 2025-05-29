def show_transformation(tokens):
    print("Starting with:", ' '.join(tokens))
    
    def can_apply_rule(i, tokens):
        if i >= len(tokens) - 1:
            return False, None
        pair = tokens[i:i+2]
        if not (pair[0].endswith('#') and pair[1].startswith('#')):
            return False, None
            
        if (pair[0] == 'A#' and pair[1] == '#A') or \
           (pair[0] == 'B#' and pair[1] == '#B'):
            return True, tokens[:i] + tokens[i+2:]
        elif pair[0] == 'A#' and pair[1] == '#B':
            return True, tokens[:i] + ['#B', 'A#'] + tokens[i+2:]
        elif pair[0] == 'B#' and pair[1] == '#A':
            return True, tokens[:i] + ['#A', 'B#'] + tokens[i+2:]
        return False, None

    step = 1
    while True:
        any_change = False
        for i in range(len(tokens)):
            can_apply, new_tokens = can_apply_rule(i, tokens)
            if can_apply:
                tokens = new_tokens
                print(f"Step {step}:", ' '.join(tokens))
                step += 1
                any_change = True
                break
        if not any_change:
            break
    
    print("Final state:", ' '.join(tokens))

initial = ['#A', 'A#', 'A#', '#A', 'A#', 'A#', 'B#', 'B#', '#A', '#B']
show_transformation(initial)