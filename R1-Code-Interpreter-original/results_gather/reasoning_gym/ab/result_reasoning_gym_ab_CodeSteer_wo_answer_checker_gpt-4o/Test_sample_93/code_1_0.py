def compute_program(tokens):
    tokens = tokens.split()
    changed = True

    while changed:
        changed = False
        new_tokens = []
        i = 0
        while i < len(tokens):
            if i < len(tokens) - 1:
                if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                    # A# #A becomes nothing
                    i += 2
                    changed = True
                elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                    # A# #B becomes #B A#
                    new_tokens.extend(['#B', 'A#'])
                    i += 2
                    changed = True
                elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                    # B# #A becomes #A B#
                    new_tokens.extend(['#A', 'B#'])
                    i += 2
                    changed = True
                elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                    # B# #B becomes nothing
                    i += 2
                    changed = True
                else:
                    new_tokens.append(tokens[i])
                    i += 1
            else:
                new_tokens.append(tokens[i])
                i += 1

        if tokens != new_tokens:
            changed = True
        tokens = new_tokens

        # Debug statement to trace the changes
        print(f"Current state: {' '.join(tokens)}")

    return ' '.join(tokens)

# Initial program
program = "B# B# #A B# #A B# #A A# A# B#"
final_state = compute_program(program)
print(f"<<<{final_state}>>>")