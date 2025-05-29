def apply_rules(program):
    tokens = program.split()
    
    def try_reduce(tokens):
        for i in range(len(tokens) - 1):
            # Rule 1: A# #A -> nothing
            if tokens[i] == "A#" and tokens[i+1] == "#A":
                return tokens[:i] + tokens[i+2:], True
            # Rule 2: B# #B -> nothing
            if tokens[i] == "B#" and tokens[i+1] == "#B":
                return tokens[:i] + tokens[i+2:], True
            # Rule 3: A# #B -> #B A#
            if tokens[i] == "A#" and tokens[i+1] == "#B":
                return tokens[:i] + ["#B", "A#"] + tokens[i+2:], True
            # Rule 4: B# #A -> #A B#
            if tokens[i] == "B#" and tokens[i+1] == "#A":
                return tokens[:i] + ["#A", "B#"] + tokens[i+2:], True
        return tokens, False

    print("Initial:", " ".join(tokens))
    steps = 0
    while steps < 100:  # Safety limit to prevent infinite loops
        tokens, changed = try_reduce(tokens)
        if not changed:
            break
        print("Step", steps + 1, ":", " ".join(tokens))
        steps += 1
    
    return " ".join(tokens)

program = "#B #A #A B# #B #A A# #A B# A#"
final = apply_rules(program)
print("\nFinal state:", final)