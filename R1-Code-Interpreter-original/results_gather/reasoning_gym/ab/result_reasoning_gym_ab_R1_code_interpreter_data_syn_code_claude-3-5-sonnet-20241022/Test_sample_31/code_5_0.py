def apply_rules(program):
    def show_state(tokens, msg=""):
        print(f"{msg}: {' '.join(tokens)}")
    
    tokens = program.split()
    show_state(tokens, "Initial")
    
    while True:
        i = 0
        changed = False
        while i < len(tokens) - 1:
            # Get current pair and check for facing ## patterns
            t1, t2 = tokens[i], tokens[i+1]
            
            # Convert #B to B# and #A to A# for matching
            pair = []
            for t in [t1, t2]:
                if t.startswith('#'):
                    pair.append(t[1:] + '#')
                else:
                    pair.append(t)
            
            # Now check the rules
            if pair[0].endswith('#') and pair[1].endswith('#'):
                if pair == ['B#', 'B#']:
                    tokens = tokens[:i] + tokens[i+2:]
                    show_state(tokens, "After B#-B# deletion")
                    changed = True
                    i = 0  # Start over from beginning
                    break
                elif pair == ['A#', 'A#']:
                    tokens = tokens[:i] + tokens[i+2:]
                    show_state(tokens, "After A#-A# deletion")
                    changed = True
                    i = 0  # Start over from beginning
                    break
            i += 1
        
        if not changed:
            # Try swapping rules only if no deletions were possible
            i = 0
            while i < len(tokens) - 1:
                t1, t2 = tokens[i], tokens[i+1]
                pair = []
                for t in [t1, t2]:
                    if t.startswith('#'):
                        pair.append(t[1:] + '#')
                    else:
                        pair.append(t)
                
                if pair[0].endswith('#') and pair[1].endswith('#'):
                    if pair == ['A#', 'B#']:
                        tokens = tokens[:i] + ['#B', 'A#'] + tokens[i+2:]
                        show_state(tokens, "After A#-B# swap")
                        changed = True
                        i = 0
                        break
                    elif pair == ['B#', 'A#']:
                        tokens = tokens[:i] + ['#A', 'B#'] + tokens[i+2:]
                        show_state(tokens, "After B#-A# swap")
                        changed = True
                        i = 0
                        break
                i += 1
        
        if not changed:
            break
    
    return tokens

program = "#B #B #A #A #B A# #B A# A# B#"
final_tokens = apply_rules(program)
print("\nFinal:", " ".join(final_tokens))