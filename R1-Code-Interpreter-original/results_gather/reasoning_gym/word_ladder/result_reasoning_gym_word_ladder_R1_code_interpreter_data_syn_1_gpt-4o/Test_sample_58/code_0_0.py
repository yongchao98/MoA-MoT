# List of valid English words (a small subset for demonstration)
valid_words = {"CLAN", "PLAN", "PLAS", "PLUS", "PUSH", "RUSH"}

# Function to check if a word is valid
def is_valid_word(word):
    return word in valid_words

# Starting and ending words
start_word = "CLAN"
end_word = "RUSH"

# Possible transformations
transformations = []

# Start with the initial word
current_word = start_word
transformations.append(current_word)

# Transformations
if is_valid_word("PLAN"):
    current_word = "PLAN"
    transformations.append(current_word)

if is_valid_word("PLAS"):
    current_word = "PLAS"
    transformations.append(current_word)

if is_valid_word("PLUS"):
    current_word = "PLUS"
    transformations.append(current_word)

if is_valid_word("PUSH"):
    current_word = "PUSH"
    transformations.append(current_word)

if is_valid_word("RUSH"):
    current_word = "RUSH"
    transformations.append(current_word)

# Print the transformations
print(",".join(transformations))