def apply_rules(program):
    def normalize_pair(t1, t2):
        # Convert #B to B# and #A to A# for comparison
        if t1.startswith('#'):
            t1 = t1[1:] + '#'
        if t2.startswith('#'):
            t2 = t2[1:] + '#'
        return t1, t2
    
    def show_state(tokens, msg=""):
        print(f"{msg}: {' '.join(tokens)}")
    
    tokens = program.split()
    show_state(tokens, "Initial")
    
    while True:
        i = 0
        changed = False
        while i < len(tokens) - 1:
            t1, t2 = normalize_pair(tokens[i], tokens[i+1])
            
            # Now check if we have a ## pattern
            if t1.endswith('#') and t2.endswith('#'):
                if (t1, t2) == ('B#', 'B#'):
                    # Delete both tokens
                    tokens = tokens[:i] + tokens[i+2:]
                    show_state(tokens, "After B#-B# deletion")
                    changed = True
                    break
                elif (t1, t2) == ('A#', 'A#'):
                    # Delete both tokens
                    tokens = tokens[:i] + tokens[i+2:]
                    show_state(tokens, "After A#-A# deletion")
                    changed = True
                    break
                elif (t1, t2) == ('A#', 'B#'):
                    # Swap and reorder
                    tokens = tokens[:i] + ['#B', 'A#'] + tokens[i+2:]
                    show_state(tokens, "After A#-B# swap")
                    changed = True
                    break
                elif (t1, t2) == ('B#', 'A#'):
                    # Swap and reorder
                    tokens = tokens[:i] + ['#A', 'B#'] + tokens[i+2:]
                    show_state(tokens, "After B#-A# swap")
                    changed = True
                    break
            i += 1
        
        if not changed:
            break
    
    return tokens

program = "#B #B #A #A #B A# #B A# A# B#"
final_tokens = apply_rules(program)
print("\nFinal:", " ".join(final_tokens))