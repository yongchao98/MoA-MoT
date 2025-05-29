def apply_rules(program):
    tokens = program.split()
    
    def try_reduce(tokens):
        # First priority: remove A# #A pairs
        for i in range(len(tokens) - 1):
            if tokens[i] == "A#" and tokens[i+1] == "#A":
                return tokens[:i] + tokens[i+2:], True
            if tokens[i] == "B#" and tokens[i+1] == "#B":
                return tokens[:i] + tokens[i+2:], True
        
        # Second priority: transform #A B# to B# #A
        for i in range(len(tokens) - 1):
            if tokens[i] == "#A" and tokens[i+1] == "B#":
                return tokens[:i] + ["B#", "#A"] + tokens[i+2:], True
            if tokens[i] == "A#" and tokens[i+1] == "#B":
                return tokens[:i] + ["#B", "A#"] + tokens[i+2:], True
            if tokens[i] == "B#" and tokens[i+1] == "#A":
                return tokens[:i] + ["#A", "B#"] + tokens[i+2:], True
            if tokens[i] == "#B" and tokens[i+1] == "A#":
                return tokens[:i] + ["A#", "#B"] + tokens[i+2:], True
        
        return tokens, False

    print("Initial:", " ".join(tokens))
    step = 0
    while step < 100:  # Safety limit
        old_tokens = tokens.copy()
        tokens, changed = try_reduce(tokens)
        if not changed:
            break
        print("Step", step + 1, ":", " ".join(tokens))
        step += 1
    
    return " ".join(tokens)

program = "#B #A #A B# #B #A A# #A B# A#"
final = apply_rules(program)
print("\nFinal state:", final)