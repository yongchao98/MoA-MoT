def apply_rules(program):
    tokens = program.split()
    
    def try_reduce(tokens):
        for i in range(len(tokens) - 1):
            # Direct rules
            if tokens[i] == "A#" and tokens[i+1] == "#A":
                return tokens[:i] + tokens[i+2:], True
            if tokens[i] == "B#" and tokens[i+1] == "#B":
                return tokens[:i] + tokens[i+2:], True
            if tokens[i] == "A#" and tokens[i+1] == "#B":
                return tokens[:i] + ["#B", "A#"] + tokens[i+2:], True
            if tokens[i] == "B#" and tokens[i+1] == "#A":
                return tokens[:i] + ["#A", "B#"] + tokens[i+2:], True
            # Pattern #A B# should transform to B# #A
            if tokens[i] == "#A" and tokens[i+1] == "B#":
                return tokens[:i] + ["B#", "#A"] + tokens[i+2:], True
            # Pattern #B A# should transform to A# #B
            if tokens[i] == "#B" and tokens[i+1] == "A#":
                return tokens[:i] + ["A#", "#B"] + tokens[i+2:], True
        return tokens, False

    print("Initial:", " ".join(tokens))
    steps = 0
    prev_states = set()
    current_state = " ".join(tokens)
    
    while current_state not in prev_states and steps < 100:
        prev_states.add(current_state)
        tokens, changed = try_reduce(tokens)
        if not changed:
            break
        current_state = " ".join(tokens)
        print("Step", steps + 1, ":", current_state)
        steps += 1
    
    return current_state

program = "#B #A #A B# #B #A A# #A B# A#"
final = apply_rules(program)
print("\nFinal state:", final)