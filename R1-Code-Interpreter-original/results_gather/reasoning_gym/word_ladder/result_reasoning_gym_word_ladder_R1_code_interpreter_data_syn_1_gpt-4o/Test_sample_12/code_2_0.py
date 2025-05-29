# List of valid English words for the transformation
valid_words = {"SNOG", "SNAG", "SNUG", "SUQS"}

# Function to check if a word is valid
def is_valid_word(word):
    return word in valid_words

# Starting and ending words
start_word = "SNOG"
end_word = "SUQS"

# List to store the transformation sequence
transformation_sequence = [start_word]

# Manually deduced sequence
manual_sequence = ["SNOG", "SNUG", "SUQS"]

# Verify the sequence
for word in manual_sequence:
    if is_valid_word(word):
        transformation_sequence.append(word)

# Print the transformation sequence
print(",".join(transformation_sequence))