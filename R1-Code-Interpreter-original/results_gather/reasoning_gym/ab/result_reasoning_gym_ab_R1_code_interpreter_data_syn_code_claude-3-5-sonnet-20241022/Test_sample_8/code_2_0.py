def can_reduce(tokens):
    for i in range(len(tokens) - 1):
        pair = tokens[i] + " " + tokens[i+1]
        if pair in ["A# #A", "B# #B", "A# #B", "B# #A"]:
            return True
    return False

def apply_single_rule(tokens, i):
    pair = tokens[i] + " " + tokens[i+1]
    if pair in ["A# #A", "B# #B"]:
        return tokens[:i] + tokens[i+2:], True
    elif pair == "A# #B":
        return tokens[:i] + ["#B", "A#"] + tokens[i+2:], True
    elif pair == "B# #A":
        return tokens[:i] + ["#A", "B#"] + tokens[i+2:], True
    return tokens, False

def apply_rules(program):
    tokens = program.split()
    step = 0
    print(f"Step {step}: {' '.join(tokens)}")
    
    while True:
        changed = False
        for i in range(len(tokens) - 1):
            new_tokens, was_changed = apply_single_rule(tokens, i)
            if was_changed:
                tokens = new_tokens
                step += 1
                print(f"Step {step}: {' '.join(tokens)}")
                changed = True
                break
        
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "#B #A #A B# B# A#"
print("Starting reduction process:")
final = apply_rules(program)
print(f"\nFinal state: {final}")