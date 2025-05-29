def apply_rules(program):
    def find_match(tokens):
        for i in range(len(tokens) - 1):
            pair = tokens[i:i+2]
            # Rule 1: A# #A -> nothing
            if pair == ['A#', '#A']:
                return i, 2, []
            # Rule 2: A# #B -> #B A#
            elif pair == ['A#', '#B']:
                return i, 2, ['#B', 'A#']
            # Rule 3: B# #A -> #A B#
            elif pair == ['B#', '#A']:
                return i, 2, ['#A', 'B#']
            # Rule 4: B# #B -> nothing
            elif pair == ['B#', '#B']:
                return i, 2, []
        return -1, 0, []

    # Convert input string to list of tokens
    tokens = program.split()
    
    while True:
        pos, length, replacement = find_match(tokens)
        if pos == -1:  # No more matches found
            break
        # Replace the matched tokens with their replacement
        tokens = tokens[:pos] + replacement + tokens[pos+length:]
    
    return ' '.join(tokens)

# Initial program
program = "#A #B #B #A A# #A #A A# A# #B"
result = apply_rules(program)
print(result)