def apply_rules(program):
    tokens = program.split()
    print("Initial:", ' '.join(tokens))
    
    step = 1
    while True:
        changed = False
        
        # First pass: do all cancellations
        i = 0
        while i < len(tokens) - 1:
            if tokens[i].endswith('A#') and tokens[i+1].startswith('#A'):
                print(f"Step {step}: Remove A# #A at positions {i},{i+1}")
                tokens.pop(i)
                tokens.pop(i)
                print("After:", ' '.join(tokens))
                changed = True
                i = 0
                step += 1
                continue
            elif tokens[i].endswith('B#') and tokens[i+1].startswith('#B'):
                print(f"Step {step}: Remove B# #B at positions {i},{i+1}")
                tokens.pop(i)
                tokens.pop(i)
                print("After:", ' '.join(tokens))
                changed = True
                i = 0
                step += 1
                continue
            i += 1
        
        # Second pass: do all swaps
        i = 0
        while i < len(tokens) - 1:
            if tokens[i].endswith('B#') and tokens[i+1].startswith('#A'):
                print(f"Step {step}: Swap B# #A at positions {i},{i+1}")
                tokens[i], tokens[i+1] = "#A", "B#"
                print("After:", ' '.join(tokens))
                changed = True
                i = 0
                step += 1
                continue
            elif tokens[i].endswith('A#') and tokens[i+1].startswith('#B'):
                print(f"Step {step}: Swap A# #B at positions {i},{i+1}")
                tokens[i], tokens[i+1] = "#B", "A#"
                print("After:", ' '.join(tokens))
                changed = True
                i = 0
                step += 1
                continue
            i += 1
            
        if not changed:
            break
    
    return ' '.join(tokens)

program = "#A #B #B B# #A A# #B B# #A #A"
final_state = apply_rules(program)
print("\nFinal state:", final_state)