def apply_rules(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == "A#" and tokens[i + 1] == "#A":
            # A# #A becomes nothing
            tokens.pop(i)
            tokens.pop(i)
        elif tokens[i] == "A#" and tokens[i + 1] == "#B":
            # A# #B becomes #B A#
            tokens[i], tokens[i + 1] = "#B", "A#"
            i += 1
        elif tokens[i] == "B#" and tokens[i + 1] == "#A":
            # B# #A becomes #A B#
            tokens[i], tokens[i + 1] = "#A", "B#"
            i += 1
        elif tokens[i] == "B#" and tokens[i + 1] == "#B":
            # B# #B becomes nothing
            tokens.pop(i)
            tokens.pop(i)
        else:
            i += 1
    return tokens

def compute_program(initial_program):
    tokens = initial_program.split()
    previous_length = -1
    while previous_length != len(tokens):
        previous_length = len(tokens)
        tokens = apply_rules(tokens)
    return tokens

initial_program = "A# #A A# B# #A #B #A #B #B A#"
final_state = compute_program(initial_program)
print(" ".join(final_state))