def apply_rules_with_steps(program):
    tokens = program.split()
    step = 1
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i+1]
            old_tokens = tokens.copy()
            
            if pair == "A# #A" or pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                print(f"Step {step}: {' '.join(old_tokens)} -> {' '.join(tokens)} (Delete rule)")
                break
            elif pair == "A# #B":
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
                print(f"Step {step}: {' '.join(old_tokens)} -> {' '.join(tokens)} (A# #B rule)")
                break
            elif pair == "B# #A":
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
                print(f"Step {step}: {' '.join(old_tokens)} -> {' '.join(tokens)} (B# #A rule)")
                break
            i += 1
        
        if not changed:
            break
        step += 1
    
    return " ".join(tokens)

program = "#A A# #B #B #A #A #B #B A# #B"
print(f"Initial: {program}")
print(f"Final: {apply_rules_with_steps(program)}")