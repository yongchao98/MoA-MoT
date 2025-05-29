def apply_rules(program):
    tokens = program.split()
    steps = [" ".join(tokens)]
    
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i+1]
            
            if pair == "A# #A":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                steps.append(" ".join(tokens))
                break
            elif pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                steps.append(" ".join(tokens))
                break
            elif pair == "A# #B":
                tokens[i] = "#B"
                tokens[i+1] = "A#"
                changed = True
                steps.append(" ".join(tokens))
                break
            elif pair == "B# #A":
                tokens[i] = "#A"
                tokens[i+1] = "B#"
                changed = True
                steps.append(" ".join(tokens))
                break
            
            i += 1
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return " ".join(tokens)

program = "A# #A #B A# #A A# A# #A B# #B"
result = apply_rules(program)
print("\nFinal result:", result)