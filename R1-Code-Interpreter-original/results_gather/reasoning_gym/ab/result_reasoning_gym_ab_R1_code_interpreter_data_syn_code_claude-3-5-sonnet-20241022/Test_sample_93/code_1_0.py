def print_steps(tokens):
    print("Initial:", " ".join(tokens))
    changed = True
    step = 1
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            t1, t2 = tokens[i], tokens[i+1]
            
            # Rule: A# #A -> nothing
            if t1 == "A#" and t2 == "#A":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                print(f"Step {step} (A# #A -> ∅):", " ".join(tokens))
                break
                
            # Rule: B# #B -> nothing
            elif t1 == "B#" and t2 == "#B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                print(f"Step {step} (B# #B -> ∅):", " ".join(tokens))
                break
                
            # Rule: A# #B -> #B A#
            elif t1 == "A#" and t2 == "#B":
                tokens[i], tokens[i+1] = tokens[i+1], tokens[i]
                changed = True
                print(f"Step {step} (A# #B -> #B A#):", " ".join(tokens))
                break
                
            # Rule: B# #A -> #A B#
            elif t1 == "B#" and t2 == "#A":
                tokens[i], tokens[i+1] = tokens[i+1], tokens[i]
                changed = True
                print(f"Step {step} (B# #A -> #A B#):", " ".join(tokens))
                break
                
            i += 1
        step += 1
    print("Final:", " ".join(tokens))

# Initial program
program = ["B#", "B#", "#A", "B#", "#A", "B#", "#A", "A#", "A#", "B#"]
print_steps(program)