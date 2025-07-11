def find_ungrammatical_sentence():
    """
    Analyzes a list of sentences in an unknown language to find the one
    that is not grammatically well-formed based on observed patterns.
    """
    sentences = [
        "Ketannet luesij gone.",
        "Ezsue esku luesij kej.",
        "Dokuje luesij ge.",
        "Kergoet dokuje otazsij ga.",
        "Dokujet esku luesij konej.",
        "Dokujet kergoe otazsij ga.",
        "Ezsuet kergoet esku otazsij kaij.",
        "Kergoet dokujet esku otazsij kosaij.",
        "Dokujet ketanne esku otazsij kaij.",
        "Ketanne esku luesij kej.",
        "Dokujet ezsuet esku otazsij kosaij.",
        "Ketannet ezsuet otazsij gosa."
    ]

    # The core grammatical rule appears in sentences with two nouns followed by the verb 'otazsij'.
    # The rule changes based on the presence of the word 'esku'.
    #
    # Deduced Rule: In sentences with two nouns followed by 'otazsij' but WITHOUT 'esku',
    # exactly one of the two nouns must have a '-t' suffix (manifesting as an '-et' ending).
    #
    # Let's check which sentence violates this rule.
    # Sentence 4: Kergoet(has -t) dokuje(no -t) -> Correct.
    # Sentence 6: Dokujet(has -t) kergoe(no -t) -> Correct.
    # Sentence 12: Ketannet(has -t) ezsuet(has -t) -> VIOLATION.

    ungrammatical_sentence_number = None

    for i, sentence in enumerate(sentences):
        words = sentence.strip(".").split()

        # Check for the pattern: two nouns followed by 'otazsij' but not 'esku'.
        if "otazsij" in words and "esku" not in words and len(words) > 2:
            noun1 = words[0]
            noun2 = words[1]

            # A simple way to check for the '-t' suffix is to see if the word ends with 'et'.
            noun1_has_t = noun1.endswith('et')
            noun2_has_t = noun2.endswith('et')

            # The rule is that exactly one should have the suffix.
            # This is equivalent to their boolean values being different (one True, one False).
            if not (noun1_has_t ^ noun2_has_t):  # Using XOR is a concise way to check this.
                ungrammatical_sentence_number = i + 1
                break  # Found the single incorrect sentence

    # The problem asks to output the number of the sentence.
    # Our logic identified sentence 12 as the one violating the rule.
    # The final output should be this number.
    if ungrammatical_sentence_number is not None:
        print(ungrammatical_sentence_number)

find_ungrammatical_sentence()
<<<12>>>