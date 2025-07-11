def solve_grammar_puzzle():
    """
    This function solves the language puzzle by deducing grammatical rules
    and finding the sentence that violates them.
    """

    # 1. Define the vocabulary and sentences
    T_FORMS = {"Ketannet", "Ezsuet", "Dokujet", "Kergoet"}
    BASE_FORMS = {"Ketanne", "Ezsue", "Dokuje", "Kergoe"}
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

    # 2. Define helper functions based on deduced rules
    def is_t_form(word):
        return word in T_FORMS

    def is_base_form(word):
        return word in BASE_FORMS

    def check_sentence_grammar(sentence):
        """
        Checks a sentence against the deduced grammatical rules.
        Returns True if grammatical, False otherwise.
        """
        words = sentence.strip('.').split()
        particle = words[-1]

        # Rule set for verb 'luesij' (1 noun argument)
        if 'luesij' in words:
            noun = words[0]
            # Rule: Particle determines noun form.
            if particle in ['gone', 'konej']:
                return is_t_form(noun)
            elif particle in ['ge', 'kej']:
                return is_base_form(noun)
            return False

        # Rule set for verb 'otazsij' (2 noun arguments)
        elif 'otazsij' in words:
            noun1 = words[0]
            noun2 = words[1]
            has_esku = 'esku' in words

            # Rule: Particle determines noun forms and 'esku' presence.
            if particle == 'ga':
                # Pattern from S4, S6: Noun-t Noun-base, no esku
                return is_t_form(noun1) and is_base_form(noun2) and not has_esku
            elif particle == 'gosa':
                # Pattern from S12: Noun-t Noun-t, no esku
                return is_t_form(noun1) and is_t_form(noun2) and not has_esku
            elif particle == 'kaij':
                # Pattern from S7 (and supported by S8, S11): Noun-t Noun-t, with esku
                return is_t_form(noun1) and is_t_form(noun2) and has_esku
            elif particle == 'kosaij':
                # Pattern from S8, S11: Noun-t Noun-t, with esku
                return is_t_form(noun1) and is_t_form(noun2) and has_esku
            return False
        
        return False

    # 3. Iterate through sentences to find the one that violates the rules
    incorrect_sentence_number = None
    for i, sentence in enumerate(sentences):
        if not check_sentence_grammar(sentence):
            # The sentence number is the index + 1
            incorrect_sentence_number = i + 1
            break
    
    # 4. Print the final answer
    if incorrect_sentence_number is not None:
        print("The number of the sentence that isn't grammatically well-formed is:")
        # The prompt asks to output each number in the final equation.
        # In this case, the result is a single number.
        print(incorrect_sentence_number)
    else:
        print("Could not determine the incorrect sentence based on the rules.")

solve_grammar_puzzle()