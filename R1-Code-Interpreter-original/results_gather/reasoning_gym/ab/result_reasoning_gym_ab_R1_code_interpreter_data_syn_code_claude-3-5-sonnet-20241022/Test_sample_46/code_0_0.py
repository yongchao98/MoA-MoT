def apply_rules(tokens):
    def find_match(tokens):
        for i in range(len(tokens)-1):
            # Rule 1: A# #A → nothing
            if tokens[i] == 'A#' and tokens[i+1] == '#A':
                return i, 2, []
            # Rule 2: A# #B → #B A#
            elif tokens[i] == 'A#' and tokens[i+1] == '#B':
                return i, 2, ['#B', 'A#']
            # Rule 3: B# #A → #A B#
            elif tokens[i] == 'B#' and tokens[i+1] == '#A':
                return i, 2, ['#A', 'B#']
            # Rule 4: B# #B → nothing
            elif tokens[i] == 'B#' and tokens[i+1] == '#B':
                return i, 2, []
        return -1, 0, []

    # Initial sequence
    sequence = tokens.split()
    
    while True:
        pos, length, replacement = find_match(sequence)
        if pos == -1:  # No more matches found
            break
            
        # Replace the matched tokens with their replacement
        sequence = sequence[:pos] + replacement + sequence[pos+length:]
        
    return ' '.join(sequence)

program = "B# A# #A #A #B A# B# #B B# #B"
result = apply_rules(program)
print(result)