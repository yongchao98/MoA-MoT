def apply_rules(tokens):
    # Returns (new_tokens, changed) where changed indicates if any rule was applied
    i = 0
    while i < len(tokens) - 1:
        # Get current pair of tokens
        pair = tokens[i] + tokens[i+1]
        
        # Check and apply rules
        if pair in ["A##A", "B##B"]:  # These pairs disappear
            tokens = tokens[:i] + tokens[i+2:]
            return tokens, True
        elif pair == "A##B":  # A# #B becomes #B A#
            tokens = tokens[:i] + ["#B", "A#"] + tokens[i+2:]
            return tokens, True
        elif pair == "B##A":  # B# #A becomes #A B#
            tokens = tokens[:i] + ["#A", "B#"] + tokens[i+2:]
            return tokens, True
        i += 1
    return tokens, False

# Initial program
program = ["#B", "B#", "A#", "A#", "#B", "#B", "#A", "A#", "#B", "A#"]

# Keep applying rules until no more changes
while True:
    program, changed = apply_rules(program)
    if not changed:
        break

# Print final state
print(" ".join(program))