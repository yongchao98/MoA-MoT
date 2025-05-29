def apply_rules(program):
    tokens = program.split()
    step = 0
    print(f"Step {step}: {' '.join(tokens)}")
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i+1]
            
            if pair == "A# #A":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                step += 1
                print(f"Step {step}: {' '.join(tokens)} (A# #A removed)")
                break
                
            elif pair == "A# #B":
                tokens[i] = "#B"
                tokens[i+1] = "A#"
                changed = True
                step += 1
                print(f"Step {step}: {' '.join(tokens)} (A# #B → #B A#)")
                break
                
            elif pair == "B# #A":
                tokens[i] = "#A"
                tokens[i+1] = "B#"
                changed = True
                step += 1
                print(f"Step {step}: {' '.join(tokens)} (B# #A → #A B#)")
                break
                
            elif pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                step += 1
                print(f"Step {step}: {' '.join(tokens)} (B# #B removed)")
                break
                
            i += 1
            
        if not changed:
            break
    
    return " ".join(tokens)

program = "A# #A A# #B #A A# A# #B #A #A"
final_state = apply_rules(program)
print("\nFinal state:", final_state)