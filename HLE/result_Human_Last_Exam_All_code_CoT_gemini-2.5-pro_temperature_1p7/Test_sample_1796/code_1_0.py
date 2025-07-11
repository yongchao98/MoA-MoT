import nltk

# This script may need to download NLTK data packs.
# We'll try to do this automatically if they are not found.
try:
    nltk.data.find('tokenizers/punkt')
    nltk.data.find('taggers/averaged_perceptron_tagger')
except LookupError:
    print("Downloading necessary NLTK data packs (punkt, averaged_perceptron_tagger)...")
    nltk.download('punkt', quiet=True)
    nltk.download('averaged_perceptron_tagger', quiet=True)
    print("Download complete.")

# The sentence to analyze
sentence = "Fat people eat accumulates."

# 1. Tokenize the sentence (split it into words)
words = nltk.word_tokenize(sentence)

# 2. Perform Part-of-Speech tagging
pos_tags = nltk.pos_tag(words)

# 3. Get the first word and its corresponding tag
first_word, first_tag = pos_tags[0]

# 4. Determine the part of speech based on the tag
part_of_speech = "Unknown"
if first_tag.startswith('JJ'):
    part_of_speech = "Adjective"
elif first_tag.startswith('NN'):
    part_of_speech = "Noun"
elif first_tag.startswith('VB'):
    part_of_speech = "Verb"
elif first_tag.startswith('RB'):
    part_of_speech = "Adverb"

# 5. Print the analysis
print(f"Sentence: \"{sentence}\"")
print(f"The first word is '{first_word}'.")
print(f"NLTK's Part-of-Speech tagger identified its role with the tag: '{first_tag}'.")
print(f"The tag '{first_tag}' represents an {part_of_speech}.")
print("\nConclusion: In this context, 'Fat' is an adjective because it modifies or describes the noun 'people'.")
