def can_apply_rule(tokens, i):
    if i >= len(tokens) - 1:
        return False
    
    t1, t2 = tokens[i], tokens[i+1]
    
    # Rule: A# #A -> nothing
    if t1 == "A#" and t2 == "#A":
        return "delete"
    
    # Rule: B# #B -> nothing
    if t1 == "B#" and t2 == "#B":
        return "delete"
    
    # Rule: A# #B -> #B A#
    if t1 == "A#" and t2 == "#B":
        return "swap"
    
    # Rule: B# #A -> #A B#
    if t1 == "B#" and t2 == "#A":
        return "swap"
    
    return False

def apply_rules(tokens):
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens):
            action = can_apply_rule(tokens, i)
            if action == "delete":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            elif action == "swap":
                tokens[i], tokens[i+1] = tokens[i+1], tokens[i]
                changed = True
                break
            i += 1
    return tokens

# Initial program
program = ["B#", "B#", "#A", "B#", "#A", "B#", "#A", "A#", "A#", "B#"]
final = apply_rules(program)
print(" ".join(final))