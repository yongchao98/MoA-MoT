# Initial sequence
sequence = ['#A', 'B#', '#B', '#B', 'A#', 'A#', 'B#', '#B', '#A', 'A#']

# Function to apply rules
def apply_rules(seq):
    i = 0
    while i < len(seq) - 1:
        if seq[i] == 'A#' and seq[i+1] == '#A':
            seq.pop(i)
            seq.pop(i)
        elif seq[i] == 'A#' and seq[i+1] == '#B':
            seq[i], seq[i+1] = '#B', 'A#'
            i += 1
        elif seq[i] == 'B#' and seq[i+1] == '#A':
            seq[i], seq[i+1] = '#A', 'B#'
            i += 1
        elif seq[i] == 'B#' and seq[i+1] == '#B':
            seq.pop(i)
            seq.pop(i)
        else:
            i += 1
    return seq

# Apply rules until no more changes
previous_sequence = []
while sequence != previous_sequence:
    previous_sequence = sequence[:]
    sequence = apply_rules(sequence)

# Print the final sequence
print(sequence)