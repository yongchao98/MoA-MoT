def compute_program(tokens):
    tokens = tokens.split()
    rewrite_occurred = True

    while rewrite_occurred:
        rewrite_occurred = False
        i = 0
        while i < len(tokens) - 1:
            if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                # A# #A becomes nothing
                tokens.pop(i)
                tokens.pop(i)
                rewrite_occurred = True
            elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                # A# #B becomes #B A#
                tokens[i], tokens[i + 1] = '#B', 'A#'
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                # B# #A becomes #A B#
                tokens[i], tokens[i + 1] = '#A', 'B#'
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                # B# #B becomes nothing
                tokens.pop(i)
                tokens.pop(i)
                rewrite_occurred = True
            else:
                i += 1

    return ' '.join(tokens)

program = "#A B# #B #B A# A# B# #B #A A#"
final_state = compute_program(program)
print(f"<<<{final_state}>>>")