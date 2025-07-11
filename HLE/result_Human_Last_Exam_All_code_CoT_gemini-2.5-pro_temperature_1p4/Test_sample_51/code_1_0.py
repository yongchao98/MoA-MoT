import textwrap

def find_inconsistent_axiom():
    """
    Analyzes axioms in the context of dependent type theory to find which one
    is inconsistent with a given structural recursion scheme.
    """

    axioms = {
        "A": "Propositional extensionality",
        "B": "Functional extensionality",
        "C": "Propositional resizing",
        "D": "Uniqueness of identity proofs",
        "E": "Proof irrelevance",
        "F": "Double-negation elimination",
        "G": "Constructive indefinite description",
        "H": "Excluded middle",
        "I": "Markov's principle"
    }

    print("--- Step 1: Understanding the Source of Inconsistency ---")
    explanation = """
    In dependent type theory, program termination is crucial for logical consistency.
    Structural recursion is a principle that guarantees termination by ensuring
    functions are only called on structurally smaller inputs. An inconsistency arises
    if an axiom allows us to bypass this guarantee.

    The vulnerability comes from IMPREDICATIVITY: the ability to define a type by
    quantifying over a collection of types that contains the one being defined.
    This allows for the creation of self-referential, paradoxical types that do not
    have a well-founded structure, which can then be used to write non-terminating
    programs. The liberal subterm relation mentioned in the prompt is known to be
    incompatible with such types.
    """
    print(textwrap.dedent(explanation))

    print("--- Step 2: Evaluating Each Axiom ---")
    # This dictionary holds the reasoning for why each axiom is or is not the culprit.
    reasoning = {
        "A": "Safe. Relates logical equivalence and equality for propositions.",
        "B": "Safe. A standard axiom about function equality.",
        "C": "DANGEROUS. This axiom is a strong form of impredicativity. It allows the universe 'Prop' to be defined in terms of itself, enabling the construction of paradoxical types (a-la Girard's paradox) that break termination.",
        "D": "Safe. Restricts the structure of equality proofs, does not add expressive power.",
        "E": "Safe. A consequence of (D) for propositions.",
        "F": "Not the cause. A classical logic principle; breaks constructivity, not termination.",
        "G": "Not the cause. A choice principle, not related to this type of paradox.",
        "H": "Not the cause. The primary classical logic axiom, same issue as (F).",
        "I": "Not the cause. A weak classical principle consistent with many systems."
    }

    for key, name in axioms.items():
        print(f"  - Axiom ({key}) {name}: {reasoning[key]}")
    
    print("\n--- Step 3: Conclusion ---")
    conclusion = """
    Propositional resizing introduces the necessary impredicativity to construct
    a non-terminating term that the specified recursion principle would accept. This
    clash between a powerful impredicative axiom and a general recursion scheme
    is a well-known source of inconsistency in type theory.
    """
    print(textwrap.dedent(conclusion))
    
    final_answer_key = "C"
    print(f"\nThe axiom inconsistent with the system is ({final_answer_key}): {axioms[final_answer_key]}.")


# Execute the analysis
find_inconsistent_axiom()
<<<C>>>