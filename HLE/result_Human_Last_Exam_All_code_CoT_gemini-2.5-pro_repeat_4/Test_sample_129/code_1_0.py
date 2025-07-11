def solve_boolean_expressions():
    """
    Calculates the number of true boolean expressions of exactly 5 symbols.

    The method is to categorize expressions by their top-level structure
    and count the number of true expressions in each category.
    """

    # Category 1: Expressions of the form (E), where E is an expression of length 3.
    # For (E) to be true, E must be a true expression of length 3.
    # From analysis, there are 5 true expressions of length 3: (T), T&T, T|T, T|F, F|T.
    # These lead to ((T)), (T&T), (T|T), (T|F), (F|T).
    cat1_paren = 5
    print(f"Category 1: Expressions of the form (E) like ((T)), (T&T), etc.")
    print(f"Count = {cat1_paren}\n")

    # Category 2: Expressions of the form !E, where E is an expression of length 4.
    # For !E to be true, E must be false. For syntactic validity, E must be parenthesized.
    # The parenthesized false expressions of length 4 are (!T) and !(T).
    # !(!T) is a valid expression of length 5 and is true.
    # !(!(T)) is length 6, so it is not counted.
    cat2_negation = 1
    print(f"Category 2: Expressions of the form !E like !(!T).")
    print(f"Count = {cat2_negation}\n")

    # Category 3: Expressions of the form A | B (lowest precedence operator is |).
    # This category is broken down by the length of A and B.
    # Case 3a: len(A)=1, len(B)=3. (e.g., T|T&T)
    #   T | B_any (6 exprs) + F | B_true (2 exprs) = 8
    cat3_a = 8
    # Case 3b: len(A)=3, len(B)=1. (e.g., T&T|T)
    #   A_any | T (6 exprs) + A_true | F (2 exprs) = 8
    cat3_b = 8
    # Case 3c: len(A)=2, len(B)=2. (e.g., !F|!T)
    #   !F|!F, !F|!T, !T|!F = 3
    cat3_c = 3
    cat3_or_total = cat3_a + cat3_b + cat3_c
    print(f"Category 3: Expressions with '|' as the main operator (A | B).")
    print(f"  - Sub-case len(A)=1, len(B)=3: {cat3_a}")
    print(f"  - Sub-case len(A)=3, len(B)=1: {cat3_b}")
    print(f"  - Sub-case len(A)=2, len(B)=2: {cat3_c}")
    print(f"Count = {cat3_a} + {cat3_b} + {cat3_c} = {cat3_or_total}\n")

    # Category 4: Expressions of the form A & B (no '|', lowest precedence operator is &).
    # Case 4a: len(A)=1, len(B)=3. (e.g., T&T&T)
    #   T & B_true = 2
    cat4_a = 2
    # Case 4b: len(A)=3, len(B)=1. (e.g., (T)&T)
    #   A_true & T = 1
    cat4_b = 1
    # Case 4c: len(A)=2, len(B)=2. (e.g., !F&!F)
    #   !F & !F = 1
    cat4_c = 1
    cat4_and_total = cat4_a + cat4_b + cat4_c
    print(f"Category 4: Expressions with '&' as the main operator (A & B).")
    print(f"  - Sub-case len(A)=1, len(B)=3: {cat4_a}")
    print(f"  - Sub-case len(A)=3, len(B)=1: {cat4_b}")
    print(f"  - Sub-case len(A)=2, len(B)=2: {cat4_c}")
    print(f"Count = {cat4_a} + {cat4_b} + {cat4_c} = {cat4_and_total}\n")

    # Final Calculation
    total_expressions = cat1_paren + cat2_negation + cat3_or_total + cat4_and_total
    print("Total Number of True Expressions:")
    print(f"{cat1_paren} (paren) + {cat2_negation} (negation) + {cat3_or_total} (or) + {cat4_and_total} (and) = {total_expressions}")

solve_boolean_expressions()