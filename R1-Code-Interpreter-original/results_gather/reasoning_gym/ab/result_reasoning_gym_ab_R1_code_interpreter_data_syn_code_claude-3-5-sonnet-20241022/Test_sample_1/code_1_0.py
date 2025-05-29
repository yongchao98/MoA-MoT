def apply_single_rule(tokens, i):
    pair = f"{tokens[i]} {tokens[i+1]}"
    if pair == "A# #A":
        return tokens[:i] + tokens[i+2:], True
    elif pair == "A# #B":
        return tokens[:i] + ["#B", "A#"] + tokens[i+2:], True
    elif pair == "B# #A":
        return tokens[:i] + ["#A", "B#"] + tokens[i+2:], True
    elif pair == "B# #B":
        return tokens[:i] + tokens[i+2:], True
    return tokens, False

def solve_with_steps(program):
    tokens = program.split()
    step = 1
    print(f"Initial: {' '.join(tokens)}")
    
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
    
    return ' '.join(tokens)

program = "B# B# #A #B B# #A #B #B #A A#"
result = solve_with_steps(program)
print(f"\nFinal result: {result}")