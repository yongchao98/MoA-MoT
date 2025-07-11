# The relevant sentences are from Noam Chomsky's 1957 book, "Syntactic Structures".
# In the section discussing grammaticality, he provides a list of examples.

# 1. Identify the last syntactically correct sentence he gives in the list.
last_correct_sentence = "the book seems interesting"
# 2. Identify the last syntactically incorrect sentence he gives.
last_incorrect_sentence = "book the seems interesting"

# 3. Extract the noun from each sentence.
# For "the book seems interesting", the noun is 'book'.
# For "book the seems interesting", the noun is also 'book'.
noun_from_correct = "book"
noun_from_incorrect = "book"

# 4. Print the results.
print(f"In Chomsky's discussion, the last syntactically correct sentence is:")
print(f"'{last_correct_sentence}'")
print(f"The noun in this sentence is: {noun_from_correct}\n")

print(f"The last syntactically incorrect sentence is:")
print(f"'{last_incorrect_sentence}'")
print(f"The noun in this sentence is: {noun_from_incorrect}")