def analyze_binding_principles():
    """
    Analyzes three sentences to determine which are ungrammatical
    due to violations of binding principles.
    """

    # --- Introduction to Binding Principles ---
    intro = (
        "Binding Theory in linguistics consists of principles that govern how different types of noun phrases can be co-referential.\n"
        "A key principle is:\n"
        "- **Principle C**: An R-expression (e.g., a name like 'Mary') must be 'free', meaning it cannot be c-commanded by a co-indexed expression.\n\n"
        "A broader view includes principles for empty categories (traces), like the **Empty Category Principle (ECP)**, which states that a trace left by movement must be properly licensed or 'bound' by its antecedent."
    )
    print(intro)
    print("="*50)

    # --- Sentence A Analysis ---
    sentence_a_analysis = (
        "Analysis of Sentence A: 'She_i likes Mary_i and Jane.'\n\n"
        "1. In this sentence, 'She' is a pronoun and 'Mary' is an R-expression.\n"
        "2. The subscript '_i' indicates that 'She' and 'Mary' are intended to refer to the same person.\n"
        "3. The pronoun 'She' c-commands the R-expression 'Mary'.\n"
        "4. **VIOLATION**: Because the R-expression 'Mary_i' is c-commanded by and co-indexed with 'She_i', it is not free. This is a direct violation of Binding Principle C.\n"
        "--> Therefore, sentence A is ungrammatical due to a binding principle violation."
    )
    print(sentence_a_analysis)
    print("="*50)

    # --- Sentence B Analysis ---
    sentence_b_analysis = (
        "Analysis of Sentence B: 'Whose does John like glasses?'\n\n"
        "1. This question is formed from a statement like 'John likes Mary's glasses' by moving the possessor 'Whose'.\n"
        "2. The resulting structure is 'Whose_i does John like [t_i glasses]', where 't_i' is the trace.\n"
        "3. This violates the 'Left Branch Condition', a constraint on movement.\n"
        "4. **VIOLATION**: This constraint is explained by the ECP. The trace 't_i' is inside a larger noun phrase and is not properly governed/licensed by its antecedent 'Whose_i', violating the ECP.\n"
        "--> Therefore, sentence B is ungrammatical due to an ECP violation, which is a principle for binding traces."
    )
    print(sentence_b_analysis)
    print("="*50)
    
    # --- Sentence C Analysis ---
    sentence_c_analysis = (
        "Analysis of Sentence C: 'Who does John like Mary and?'\n\n"
        "1. This question is formed from 'John likes Mary and himself' by questioning 'himself' and moving it.\n"
        "2. The resulting structure is 'Who_i does John like [Mary and t_i]', where 't_i' is the trace.\n"
        "3. This violates the 'Coordinate Structure Constraint', which forbids moving an element out of a coordinated phrase.\n"
        "4. **VIOLATION**: This constraint is also explained by the ECP. The trace 't_i' is inside a coordinate structure and is not properly governed/licensed by its antecedent 'Who_i', violating the ECP.\n"
        "--> Therefore, sentence C is ungrammatical due to an ECP violation."
    )
    print(sentence_c_analysis)
    print("="*50)

    # --- Conclusion ---
    conclusion = (
        "CONCLUSION:\n"
        "- Sentence A violates Binding Principle C.\n"
        "- Sentences B and C violate the Empty Category Principle (ECP), which is part of the binding system for traces.\n"
        "- Thus, all three sentences (A, B, and C) are ungrammatical due to violations of binding principles in their broad sense."
    )
    print(conclusion)

if __name__ == '__main__':
    analyze_binding_principles()