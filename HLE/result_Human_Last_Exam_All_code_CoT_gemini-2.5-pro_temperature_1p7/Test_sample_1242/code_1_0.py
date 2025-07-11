# This script identifies and displays the group of linguistic features
# known from research to be strong predictors of word complexity.

# Based on linguistic and psycholinguistic studies, the most contributive
# features for predicting word complexity are a mix of structural and
# psycholinguistic properties of a word.

# The identified group of features is:
# 1. Word Length: A basic structural measure of complexity.
# 2. Word Familiarity Rating: A psycholinguistic measure of how common a word is.
# 3. Syllable Count: A more refined structural measure than word length.
# 4. Concreteness Rating: A psycholinguistic measure of how abstract a word's meaning is.

correct_features = [
    "word length",
    "word familiarity rating",
    "syllable count",
    "concreteness rating"
]

print("The linguistic features strongly correlated with word complexity are:")
for feature in correct_features:
    print(f"- {feature}")
