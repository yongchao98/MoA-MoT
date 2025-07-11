def find_max_fixed_points():
    """
    This function explains the reasoning to find the maximum number of fixed points
    for the composition of two degree-3 polynomials with positive derivatives.
    """

    # Introduction to the problem
    print("Problem: Let f and g be polynomials of degree 3 such that f'(x) > 0 and g'(x) > 0 for all x.")
    print("What is the maximum number of fixed points that f(g(x)) can have?")
    print("-" * 70)

    # Step 1: Formulating the fixed-point equation
    print("Step 1: The fixed-point equation")
    print("A fixed point 'x' of a function h(x) is a value such that h(x) = x.")
    print("For our problem, h(x) = f(g(x)). So we need to find the number of solutions to:")
    print("f(g(x)) = x")
    print("This can be rewritten as finding the roots of the equation: f(g(x)) - x = 0.")
    print("\n")

    # Step 2: Determining the degree of the polynomial
    print("Step 2: Degree of the fixed-point polynomial")
    deg_f = 3
    deg_g = 3
    deg_h = deg_f * deg_g
    print(f"f(x) is a polynomial of degree {deg_f}.")
    print(f"g(x) is a polynomial of degree {deg_g}.")
    print(f"The composition h(x) = f(g(x)) is a polynomial of degree {deg_f} * {deg_g} = {deg_h}.")
    print(f"Therefore, the fixed-point equation f(g(x)) - x = 0 is a polynomial equation of degree {deg_h}.")
    print("\n")

    # Step 3: Finding the upper bound for the number of roots
    print("Step 3: Upper bound on the number of fixed points")
    print(f"By the Fundamental Theorem of Algebra, a polynomial of degree {deg_h} can have at most {deg_h} real roots.")
    print(f"This means the maximum possible number of fixed points is {deg_h}.")
    print("\n")

    # Step 4: Arguing that the upper bound is achievable
    print("Step 4: Showing that 9 fixed points are achievable")
    print("Let k(x) = f(g(x)) - x. For this 9th-degree polynomial to have 9 distinct real roots, its derivative, k'(x), must have 8 distinct real roots.")
    print("Let's analyze k'(x):")
    print("k'(x) = d/dx [f(g(x)) - x] = f'(g(x)) * g'(x) - 1")
    print("f'(x) is a degree 2 polynomial, and g(x) is a degree 3 polynomial.")
    print(" - The degree of f'(g(x)) is 2 * 3 = 6.")
    print(" - The degree of g'(x) is 2.")
    print(" - The degree of their product, f'(g(x))g'(x), is 6 + 2 = 8.")
    print("So, k'(x) is an 8th-degree polynomial. A polynomial of degree 8 can have up to 8 real roots.")
    print("\nThe conditions f'(x)>0 and g'(x)>0 are achievable (e.g., by ensuring the quadratic derivatives have no real roots and open upwards).")
    print("It is possible to choose the coefficients of f and g such that:")
    print("  1. The derivative conditions f' > 0 and g' > 0 are met.")
    print("  2. The 8th-degree polynomial k'(x) has 8 distinct real roots. This creates 8 local extrema for k(x).")
    print("  3. By adjusting the constant term of f(x) (which does not affect f' or g'), we can shift k(x) vertically.")
    print("We can position the graph of k(x) so that its 4 local maxima are positive and its 4 local minima are negative.")
    print("This configuration guarantees 9 distinct real roots for k(x) = 0 by the Intermediate Value Theorem.")
    print("\n")

    # Final Conclusion
    max_points = 9
    print("Conclusion:")
    print(f"The maximum number of fixed points is the maximum number of real roots for the 9th-degree polynomial equation f(g(x)) - x = 0.")
    print("As argued, it is possible to construct such functions f and g for which 9 distinct roots exist.")
    print(f"The maximum number is {max_points}.")

if __name__ == '__main__':
    find_max_fixed_points()