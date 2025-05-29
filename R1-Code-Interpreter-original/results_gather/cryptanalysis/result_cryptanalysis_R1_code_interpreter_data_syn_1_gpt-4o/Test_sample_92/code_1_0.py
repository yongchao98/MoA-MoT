# Remaining possible letters after deducing X is correct but in the wrong position
possible_letters = {'R', 'Z', 'E', 'N', 'Q', 'Y'}

# Since X is correct but in the wrong position, it should be the last letter
# We need to find the correct letter for the third position

# Test each remaining letter for the third position
for letter in possible_letters:
    print(f"Testing combination: 2, 6, {letter}, X")