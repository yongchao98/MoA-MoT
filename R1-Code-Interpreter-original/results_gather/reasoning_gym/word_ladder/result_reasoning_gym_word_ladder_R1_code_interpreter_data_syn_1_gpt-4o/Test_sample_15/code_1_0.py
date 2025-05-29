# Manually constructed sequence
sequence = ["BRAD", "BEAD", "BEAM", "BOAM", "POAM", "POEM"]

# Check if each word in the sequence is valid
word_set = set(["BRAD", "BEAD", "BEAM", "BOAM", "POAM", "POEM"])

# Verify the sequence
is_valid_sequence = all(word in word_set for word in sequence)

# Print the sequence if valid
if is_valid_sequence:
    print(",".join(sequence))
else:
    print("No valid sequence found.")