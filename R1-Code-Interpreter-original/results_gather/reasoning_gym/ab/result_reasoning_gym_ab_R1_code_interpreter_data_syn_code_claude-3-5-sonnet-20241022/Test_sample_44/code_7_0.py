def apply_rules(program):
    tokens = program.split()
    print("Initial:", ' '.join(tokens))
    
    step = 1
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # First priority: handle cancellations
            if tokens[i] == "A#" and tokens[i+1] == "#A":
                print(f"Step {step}: Cancel A# #A at {i},{i+1}")
                tokens.pop(i)
                tokens.pop(i)
                print("After:", ' '.join(tokens))
                changed = True
                i = 0
                step += 1
                continue
            
            if tokens[i] == "B#" and tokens[i+1] == "#B":
                print(f"Step {step}: Cancel B# #B at {i},{i+1}")
                tokens.pop(i)
                tokens.pop(i)
                print("After:", ' '.join(tokens))
                changed = True
                i = 0
                step += 1
                continue
            
            # Second priority: handle transformations
            if tokens[i] == "B#" and tokens[i+1] == "#A":
                print(f"Step {step}: Transform B# #A to #A B# at {i},{i+1}")
                tokens[i] = "#A"
                tokens[i+1] = "B#"
                print("After:", ' '.join(tokens))
                changed = True
                i = 0
                step += 1
                continue
            
            if tokens[i] == "A#" and tokens[i+1] == "#B":
                print(f"Step {step}: Transform A# #B to #B A# at {i},{i+1}")
                tokens[i] = "#B"
                tokens[i+1] = "A#"
                print("After:", ' '.join(tokens))
                changed = True
                i = 0
                step += 1
                continue
            
            i += 1
        
        # Extra check: look for B# #A pattern in the entire sequence
        for i in range(len(tokens)):
            if tokens[i] == "B#":
                for j in range(i+1, len(tokens)):
                    if tokens[j] == "#A":
                        print(f"Step {step}: Found distant B# #A at {i},{j}")
                        # Perform necessary swaps to bring them together
                        changed = True
                        break
                if changed:
                    break
        
        if not changed:
            break
    
    return ' '.join(tokens)

program = "#A #B #B B# #A A# #B B# #A #A"
final_state = apply_rules(program)
print("\nFinal state:", final_state)