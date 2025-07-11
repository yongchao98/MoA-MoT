def solve_type_theory_problem():
    """
    Analyzes axioms in dependent type theory to find which is inconsistent
    with structural recursion.
    """

    axioms = {
        'A': {'name': 'Propositional extensionality',
              'description': 'If P and Q are logically equivalent, they are equal as types. Generally considered safe.',
              'is_inconsistent': False},
        'B': {'name': 'Functional extensionality',
              'description': 'If two functions behave identically, they are equal. Generally considered safe.',
              'is_inconsistent': False},
        'C': {'name': 'Propositional resizing',
              'description': 'A proposition in a high universe can be "resized" to an equivalent one in a low universe.',
              'is_inconsistent': True},
        'D': {'name': 'Uniqueness of identity proofs',
              'description': 'Any two proofs of the same equality are equal. Consistent in many systems.',
              'is_inconsistent': False},
        'E': {'name': 'Proof irrelevance',
              'description': 'Any two proofs of the same proposition are equal. Consistent.',
              'is_inconsistent': False},
        'F': {'name': 'Double-negation elimination',
              'description': 'A principle of classical logic. Not the source of this kind of paradox.',
              'is_inconsistent': False},
        'G': {'name': 'Constructive indefinite description',
              'description': 'A form of the Axiom of Choice. Generally consistent in constructive forms.',
              'is_inconsistent': False},
        'H': {'name': 'Excluded middle',
              'description': 'A strong principle of classical logic. Not the source of this paradox.',
              'is_inconsistent': False},
        'I': {'name': 'Markov\'s principle',
              'description': 'A weak classical principle. Not the source of this paradox.',
              'is_inconsistent': False}
    }

    print("Analyzing which axiom is inconsistent with structural recursion in dependent type theory...")
    print("="*70)
    print("Context: Structural recursion guarantees that programs terminate. An axiom that allows creating non-terminating programs breaks the logical consistency of the system.")
    print("The primary cause of such paradoxes in type theory is 'impredicativity', where a type is defined by quantifying over a collection that includes the new type itself.")

    culprit_key = None
    for key, info in axioms.items():
        if info['is_inconsistent']:
            culprit_key = key
            break

    if not culprit_key:
        print("\nAnalysis complete: No inconsistent axiom found in the list.")
        return

    culprit_info = axioms[culprit_key]
    print(f"\nIdentified Inconsistent Axiom: [{culprit_key}] {culprit_info['name']}")
    print("-" * 70)
    print("Explanation:\n")
    print("Modern type theories avoid Girard's paradox by using a hierarchy of universes (Type₀, Type₁, ...), where definitions quantifying over Typeₙ result in a type in a strictly larger universe, Typeₙ₊₁ This is called 'predicativity'.")
    print("\nThe axiom of 'Propositional Resizing' undermines this safety mechanism. It states that for any proposition P living in a large universe (e.g., Typeₙ), there exists an equivalent proposition P' in the smallest universe of propositions (Type₀).")
    print("\nThis allows one to 'squash' large, complex propositions into the base universe, effectively reintroducing impredicativity. With this power, one can construct a type W (a type of 'ordinals' or 'well-founded trees') by quantifying over all propositions in Type₀, and then, using resizing, prove that a statement *about* W is itself a proposition in Type₀. This allows constructing a member of W from W itself, leading to a contradiction where an 'ordinal' is smaller than itself. This violates the well-foundedness required for structural recursion.")

    print("\nConclusion:")
    print("The other axioms on the list are known to be consistent with structural recursion in predicative type theory. The classical logic axioms are a separate concern from the structural paradoxes caused by impredicativity. Therefore, Propositional Resizing is the axiom that is inconsistent with the given setup.")


solve_type_theory_problem()