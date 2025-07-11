# Noam Chomsky's "Syntactic Structures" (1957) provides sentences to
# demonstrate the independence of syntax from semantics.

# The last syntactically correct sentence given in his list is:
last_correct_sentence = "have you a book on modern music?"

# The last syntactically incorrect sentence given is:
last_incorrect_sentence = "read you a book on modern music?"

# We need to find the two nouns used in these sentences.
# In "have you a book on modern music?", the nouns are 'book' and 'music'.
# In "read you a book on modern music?", the nouns are also 'book' and 'music'.

noun1 = "book"
noun2 = "music"

print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
print(f"The last syntactically incorrect sentence is: '{last_incorrect_sentence}'")
print(f"The two nouns used in these sentences are: '{noun1}' and '{noun2}'.")