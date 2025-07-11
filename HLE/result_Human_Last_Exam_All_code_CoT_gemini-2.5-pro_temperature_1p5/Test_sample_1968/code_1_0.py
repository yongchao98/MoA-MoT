def prove_non_existence_of_function():
    """
    This script outlines the proof that for any infinite cardinal κ,
    there is no function f: [κ⁺]² → κ such that for every subset x ⊆ κ⁺
    of order type κ+1, the image f''[x]² has cardinality κ.

    The proof is by contradiction.
    """

    print("--- The Proof ---")
    print("\nStep 1: State the assumption for contradiction.")
    print("Assume such a function f exists for an infinite cardinal κ. This means:")
    print("For EVERY subset x ⊆ κ⁺ with order type κ+1, we have |f''[x]²| = κ.")

    print("\nStep 2: Apply the Canonization Lemma.")
    print("The Erdős-Hajnal-Rado Canonization Lemma applies to any function on pairs from a regular cardinal.")
    print("Since κ⁺ is the successor of a cardinal, it is always a regular cardinal.")
    print("The lemma guarantees there is a large subset A ⊆ κ⁺ (with |A| = κ⁺) on which f has a simple 'canonical' form.")

    print("\nStep 3: Analyze the canonical forms on A.")
    print("For α < β, both in A, f({α, β}) must be one of:")
    print("  1. Constant: f({α, β}) = c")
    print("  2. Projection on minimum: f({α, β}) = g(α)")
    print("  3. Projection on maximum: f({α, β}) = g(β)")
    print("  4. Injective")
    print("Case 4 is impossible because the domain [A]² has size κ⁺, but the codomain κ has size κ, and κ⁺ > κ.")

    print("\nStep 4: Reduce the remaining cases to a constant function.")
    print("In Cases 2 and 3, g is a function from A (size κ⁺) to κ (size κ).")
    print("By the pigeonhole principle, there must be a value γ ∈ κ for which the preimage g⁻¹({γ}) is large.")
    print("This means there's a subset A' ⊆ A of size κ⁺ where g is constant.")
    print("For any α, β ∈ A' (with α < β), f({α, β}) will equal this constant value γ.")
    print("So, in all possible cases (1, 2, and 3), we can find a set A' ⊆ κ⁺ of size κ⁺ on which f is constant.")

    print("\nStep 5: Construct the contradictory set x.")
    print("Since |A'| = κ⁺, and κ+1 ≤ κ⁺, we can pick a subset x ⊆ A' that has an order type of κ+1.")

    print("\nStep 6: Derive the contradiction.")
    print("  (a) By our initial assumption on f, because x has order type κ+1, the image size must be κ: |f''[x]²| = κ.")
    print("  (b) By our construction, x is a subset of A', where f is constant. Therefore, the image of f on pairs from x is a single element. The image size is 1: |f''[x]²| = 1.")

    print("\nEquating the two results leads to a contradiction:")
    
    kappa_symbol = "κ"
    derived_value = 1
    
    print("Final Equation: {} = {}".format(derived_value, kappa_symbol))
    
    print("\nThis is a contradiction, as κ is an infinite cardinal and cannot equal 1.")

    print("\n--- Conclusion ---")
    print("The initial assumption must be false. Such a function can never exist, regardless of the choice of the infinite cardinal κ.")

prove_non_existence_of_function()