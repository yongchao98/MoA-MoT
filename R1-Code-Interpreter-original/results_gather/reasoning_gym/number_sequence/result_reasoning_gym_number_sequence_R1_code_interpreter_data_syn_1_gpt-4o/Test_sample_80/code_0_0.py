# Define the sequence
sequence = [1, 0, 0, 0, 0]

# Check the pattern
def predict_next(sequence):
    # If the first element is 1 and all others are 0, predict the next as 0
    if sequence[0] == 1 and all(x == 0 for x in sequence[1:]):
        return 0
    else:
        return None

# Predict the next term
next_term = predict_next(sequence)
print(next_term)