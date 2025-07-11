def solve_curve_reduction():
    """
    Analyzes the stable reduction of the given curve at p=2.
    """
    # Coefficients of y^2 = c5*x^5 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x
    # from y^2=8*x+1*x^2+4*x^3+4*x^4+8*x^5
    c = {
        5: 8,
        4: 4,
        3: 4,
        2: 1,
        1: 8
    }

    print("Original curve equation: y^2 = 8*x^5 + 4*x^4 + 4*x^3 + 1*x^2 + 8*x")
    print("Reducing this model modulo 2 gives y^2 = x^2, which is singular.\n")
    
    # We found the transformation n=3 works.
    # x = 2^3 * x_new, y = 2^3 * y_new
    n = 3
    p = 2
    
    print(f"Applying the transformation x = {p**n}*x' and y = {p**n}*y':")

    # Calculate new coefficients for y_new^2 = c'_5*x_new^5 + ...
    # c'_i = c_i * 2^(i*n - 2*n)
    c_new = {i: c[i] * (p ** (i * n - 2 * n)) for i in c}

    print("The transformed equation is:")
    equation_str = "y'^2 = "
    terms = []
    # Print the equation with its coefficients
    for i in sorted(c_new.keys(), reverse=True):
        terms.append(f"{c_new[i]}*x'^{i}")
    equation_str += " + ".join(terms)
    print(equation_str)
    
    print("\nThe coefficients of the new equation are:")
    for i in sorted(c_new.keys(), reverse=True):
        print(f"c'_{i}: {c_new[i]}")

    print("\nReducing the new equation modulo 2:")
    reduced_eq = "y'^2 = x'^2 + x'"
    print(reduced_eq)

    print("\nTo check for singularities, we check the partial derivatives of G(x', y') = y'^2 - x'^2 - x'.")
    print("In F_2, G(x', y') = y'^2 + x'^2 + x'.")
    print("dG/dx' = 2x' + 1 = 1 (mod 2)")
    print("dG/dy' = 2y' = 0 (mod 2)")
    print("Since dG/dx' is never zero, there are no singular points.")
    
    print("\nThe reduced curve is smooth. This means the original curve has good reduction at p=2.")
    print("The stable reduction is this smooth curve, which has no double points.")
    
    num_double_points = 0
    print(f"\nNumber of double points: {num_double_points}")


solve_curve_reduction()