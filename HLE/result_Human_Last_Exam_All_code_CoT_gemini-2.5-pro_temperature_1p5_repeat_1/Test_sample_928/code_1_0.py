# In Noam Chomsky's "Syntactic Structures", the section discussing "Colorless green ideas sleep furiously"
# ends with a pair of sentences to demonstrate the difference between syntactically correct and incorrect sentences.
# Last correct sentence: "I saw a fragile whale"
# Last incorrect sentence: "I saw a fragile of"

# We identify the nouns used in these two sentences. Pronouns are considered a class of noun.
# Nouns in the correct sentence are "I" and "whale".
# The noun in the incorrect sentence is "I".
# The two unique nouns used across both sentences are "I" and "whale".

final_nouns = ["I", "whale"]

print("The two nouns used in the last syntactically correct and incorrect sentences are:")
for noun in final_nouns:
  print(noun)