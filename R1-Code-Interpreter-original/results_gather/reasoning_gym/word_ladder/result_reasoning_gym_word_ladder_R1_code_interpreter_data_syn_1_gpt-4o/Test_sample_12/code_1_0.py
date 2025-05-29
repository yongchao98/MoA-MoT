# List of valid English words for the transformation
valid_words = {"SNOG", "SNAG", "SLAG", "SLUG", "SLOG", "SLOP", "SUQS", "SUQS"}

# Function to check if a word is valid
def is_valid_word(word):
    return word in valid_words

# Starting and ending words
start_word = "SNOG"
end_word = "SUQS"

# List to store the transformation sequence
transformation_sequence = [start_word]

# Current word to transform
current_word = start_word

# Transform the word
for i in range(len(current_word)):
    for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        # Create a new word by changing one letter
        new_word = current_word[:i] + letter + current_word[i+1:]
        # Check if the new word is valid and not already in the sequence
        if is_valid_word(new_word) and new_word not in transformation_sequence:
            transformation_sequence.append(new_word)
            current_word = new_word
            if current_word == end_word:
                break
    if current_word == end_word:
        break

# Print the transformation sequence
print(",".join(transformation_sequence))