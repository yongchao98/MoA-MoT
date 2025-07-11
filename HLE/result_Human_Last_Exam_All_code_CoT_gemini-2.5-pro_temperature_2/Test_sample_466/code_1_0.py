import math

def count_adjunctions():
    """
    Calculates the number of internal adjunctions based on the derived formula.
    The number of adjunctions is equal to the number of commuting pairs of automorphisms
    in GL(3, F_11), which is |GL(3, F_11)| * k(GL(3, F_11)).
    """
    q = 11
    n = 3

    # Step 1: Calculate the size of GL(3, F_11)
    # |GL(n,q)| = (q^n - 1)(q^n - q)...(q^n - q^(n-1))
    print("Calculating the size of GL(3, F_11):")
    term1 = q**n - q**0
    term2 = q**n - q**1
    term3 = q**n - q**2
    print(f"  Term 1 (q^n - 1): {term1}")
    print(f"  Term 2 (q^n - q): {term2}")
    print(f"  Term 3 (q^n - q^2): {term3}")
    
    gl_size = term1 * term2 * term3
    print(f"Size of GL(3, {q}) = {term1} * {term2} * {term3} = {gl_size}\n")

    # Step 2: Calculate the number of conjugacy classes k(GL(3, F_11))
    # Number of non-zero elements in F_11 is q-1 = 10.
    
    print("Calculating the number of conjugacy classes k(GL(3, F_11)):")
    # Case 1: 3 distinct eigenvalues from F_11^*
    # This corresponds to choosing 3 distinct values from 10.
    c1 = math.comb(q - 1, 3)
    print(f"  Classes from 3 distinct eigenvalues: {c1}")
    
    # Case 2: 2 distinct eigenvalues (one with multiplicity 2)
    # Choose the two eigenvalues in order, times 2 partitions for the dim-2 eigenspace.
    c2 = (q - 1) * (q - 2) * 2
    print(f"  Classes from 2 distinct eigenvalues: {c2}")

    # Case 3: 1 distinct eigenvalue (with multiplicity 3)
    # Choose the eigenvalue, times 3 partitions of 3.
    c3 = (q - 1) * 3
    print(f"  Classes from 1 distinct eigenvalue: {c3}")

    # Case 4: 1 linear factor and 1 irreducible quadratic factor.
    # Number of monic irreducible polys of degree 2 over F_q is (q^2-q)/2
    num_irred_quad = (q**2 - q) // 2
    c4 = (q - 1) * num_irred_quad
    print(f"  Classes from 1 linear + 1 irred. quadratic factor: {c4}")

    # Case 5: 1 irreducible cubic factor.
    # Number of monic irreducible polys of degree 3 over F_q is (q^3-q)/3
    num_irred_cubic = (q**3 - q) // 3
    c5 = num_irred_cubic
    print(f"  Classes from 1 irred. cubic factor: {c5}")
    
    num_classes = c1 + c2 + c3 + c4 + c5
    print(f"Total number of conjugacy classes k(GL(3, {q})) = {c1} + {c2} + {c3} + {c4} + {c5} = {num_classes}\n")
    
    # Final Calculation
    total_adjunctions = gl_size * num_classes
    print("Final Calculation:")
    print(f"  Number of adjunctions = |GL(3, {q})| * k(GL(3, {q}))")
    print(f"  = {gl_size} * {num_classes}")
    print(f"  = {total_adjunctions}")

count_adjunctions()