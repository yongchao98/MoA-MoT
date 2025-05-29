def can_apply_rule(token1, token2):
    # Check if two tokens can be rewritten
    if token1 == "A#" and token2 == "#A": return True, []
    if token1 == "A#" and token2 == "#B": return True, ["#B", "A#"]
    if token1 == "B#" and token2 == "#A": return True, ["#A", "B#"]
    if token1 == "B#" and token2 == "#B": return True, []
    return False, None

def apply_rules_once(tokens):
    # Try to apply one rule
    for i in range(len(tokens) - 1):
        can_apply, result = can_apply_rule(tokens[i], tokens[i+1])
        if can_apply:
            return tokens[:i] + result + tokens[i+2:], True
    return tokens, False

def solve_program(initial):
    tokens = initial.split()
    print("Initial:", " ".join(tokens))
    
    step = 1
    while True:
        tokens, changed = apply_rules_once(tokens)
        if changed:
            print(f"Step {step}:", " ".join(tokens))
            step += 1
        else:
            break
    
    return " ".join(tokens)

initial = "#B B# #B #B #A A# #B #A A# A#"
final = solve_program(initial)
print("\nFinal state:", final)