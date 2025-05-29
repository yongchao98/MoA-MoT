def apply_rules(program):
    tokens = program.split()
    i = 0
    print("Initial:", " ".join(tokens))
    step = 1
    
    while i < len(tokens) - 1:
        pair = f"{tokens[i]} {tokens[i+1]}"
        original_tokens = tokens.copy()
        
        if pair == "A# #A":
            tokens.pop(i)
            tokens.pop(i)
            i = max(0, i-1)
            print(f"Step {step} (A# #A -> ∅):", " ".join(tokens))
        elif pair == "A# #B":
            tokens[i] = "#B"
            tokens[i+1] = "A#"
            print(f"Step {step} (A# #B -> #B A#):", " ".join(tokens))
            continue
        elif pair == "B# #A":
            tokens[i] = "#A"
            tokens[i+1] = "B#"
            print(f"Step {step} (B# #A -> #A B#):", " ".join(tokens))
            continue
        elif pair == "B# #B":
            tokens.pop(i)
            tokens.pop(i)
            i = max(0, i-1)
            print(f"Step {step} (B# #B -> ∅):", " ".join(tokens))
        else:
            i += 1
        
        if tokens != original_tokens:
            step += 1
    
    return " ".join(tokens)

program = "B# B# #B A# #A B# B# #A #A #A"
result = apply_rules(program)
print("\nFinal result:", result)