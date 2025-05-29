def apply_rules(program):
    tokens = program.split()
    print("Initial:", ' '.join(tokens))
    
    step = 1
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            current = tokens[i]
            next_token = tokens[i+1]
            
            # Check for B# #A pattern specifically
            if current == "B#" and next_token == "#A":
                print(f"Step {step}: Found B# #A at positions {i},{i+1}")
                tokens[i] = "#A"
                tokens[i+1] = "B#"
                print("After:", ' '.join(tokens))
                changed = True
                i = 0  # Start over from beginning
                step += 1
                continue
            
            # Check for A# #A pattern
            if current == "A#" and next_token == "#A":
                print(f"Step {step}: Remove A# #A at positions {i},{i+1}")
                tokens.pop(i)
                tokens.pop(i)
                print("After:", ' '.join(tokens))
                changed = True
                i = 0
                step += 1
                continue
            
            # Check for B# #B pattern
            if current == "B#" and next_token == "#B":
                print(f"Step {step}: Remove B# #B at positions {i},{i+1}")
                tokens.pop(i)
                tokens.pop(i)
                print("After:", ' '.join(tokens))
                changed = True
                i = 0
                step += 1
                continue
            
            # Check for A# #B pattern
            if current == "A#" and next_token == "#B":
                print(f"Step {step}: Swap A# #B at positions {i},{i+1}")
                tokens[i] = "#B"
                tokens[i+1] = "A#"
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