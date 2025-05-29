def show_transformation(tokens):
    # Helper function to show each step
    print(' '.join(tokens))

def apply_rules(program):
    tokens = program.split()
    show_transformation(tokens)
    
    step = 1
    while True:
        changed = False
        for i in range(len(tokens) - 1):
            # First, check for cancellations (they take priority)
            if tokens[i].endswith('A#') and tokens[i+1].startswith('#A'):
                print(f"\nStep {step}: Remove A# #A at position {i}")
                tokens.pop(i)
                tokens.pop(i)
                show_transformation(tokens)
                changed = True
                break
            elif tokens[i].endswith('B#') and tokens[i+1].startswith('#B'):
                print(f"\nStep {step}: Remove B# #B at position {i}")
                tokens.pop(i)
                tokens.pop(i)
                show_transformation(tokens)
                changed = True
                break
            # Then check for swaps
            elif tokens[i].endswith('A#') and tokens[i+1].startswith('#B'):
                print(f"\nStep {step}: Swap A# #B at position {i}")
                tokens[i], tokens[i+1] = "#B", "A#"
                show_transformation(tokens)
                changed = True
                break
            elif tokens[i].endswith('B#') and tokens[i+1].startswith('#A'):
                print(f"\nStep {step}: Swap B# #A at position {i}")
                tokens[i], tokens[i+1] = "#A", "B#"
                show_transformation(tokens)
                changed = True
                break
        
        if not changed:
            break
        step += 1
    
    return ' '.join(tokens)

program = "#A #B #B B# #A A# #B B# #A #A"
final_state = apply_rules(program)
print("\nFinal state:", final_state)