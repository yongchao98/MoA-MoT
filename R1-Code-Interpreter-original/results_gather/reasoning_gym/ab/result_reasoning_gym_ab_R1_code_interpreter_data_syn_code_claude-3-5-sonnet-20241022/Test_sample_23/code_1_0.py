def apply_single_rule(tokens, i):
    pair = tokens[i] + " " + tokens[i+1]
    
    if pair == "A# #A":  # Rule 1
        return tokens[:i] + tokens[i+2:], True
    elif pair == "A# #B":  # Rule 2
        return tokens[:i] + ["#B", "A#"] + tokens[i+2:], True
    elif pair == "B# #A":  # Rule 3
        return tokens[:i] + ["#A", "B#"] + tokens[i+2:], True
    elif pair == "B# #B":  # Rule 4
        return tokens[:i] + tokens[i+2:], True
    
    return tokens, False

def apply_rules_with_steps(program):
    tokens = program.split()
    step = 1
    
    while True:
        changed = False
        for i in range(len(tokens) - 1):
            new_tokens, was_changed = apply_single_rule(tokens, i)
            if was_changed:
                tokens = new_tokens
                print(f"Step {step}: {' '.join(tokens)}")
                step += 1
                changed = True
                break
                
        if not changed:
            break
    
    return " ".join(tokens)

program = "A# B# B# #A A# B# #A A# #A #B"
print("Initial:", program)
result = apply_rules_with_steps(program)
print("\nFinal:", result)