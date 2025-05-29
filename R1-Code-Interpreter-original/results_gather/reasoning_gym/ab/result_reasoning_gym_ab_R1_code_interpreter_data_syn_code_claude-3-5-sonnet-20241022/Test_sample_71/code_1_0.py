def apply_rules(tokens):
    i = 0
    while i < len(tokens) - 1:
        # Check if we have a facing ## pattern
        current = tokens[i]
        next_token = tokens[i+1]
        
        # Extract the # position (start or end) for each token
        curr_hash = 0 if current.startswith('#') else 1
        next_hash = 0 if next_token.startswith('#') else 1
        
        # If # symbols are adjacent (one at end, one at start)
        if curr_hash == 1 and next_hash == 0:
            pair = current + next_token
            if pair in ["A##A", "B##B"]:  # Disappear rule
                tokens = tokens[:i] + tokens[i+2:]
                return tokens, True
            elif pair == "A##B":  # A# #B -> #B A#
                tokens = tokens[:i] + ["#B", "A#"] + tokens[i+2:]
                return tokens, True
            elif pair == "B##A":  # B# #A -> #A B#
                tokens = tokens[:i] + ["#A", "B#"] + tokens[i+2:]
                return tokens, True
        i += 1
    return tokens, False

# Initial program
program = ["#B", "B#", "A#", "A#", "#B", "#B", "#A", "A#", "#B", "A#"]

# Keep applying rules until no more changes
steps = 0
while True:
    program, changed = apply_rules(program)
    if not changed:
        break
    steps += 1

print(" ".join(program))