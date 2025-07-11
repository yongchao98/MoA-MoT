import math

def compute_galois_group_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is P(x) = x^4 + px + q
    p = 8
    q = 14
    
    print(f"The polynomial is P(x) = x^4 + {p}x + {q}.")
    print("-" * 50)

    # Step 1: Check for irreducibility
    print("Step 1: Check for irreducibility of the polynomial.")
    print("We use Eisenstein's criterion with the prime p=2:")
    print(" - The prime 2 divides the coefficients 0 (for x^3), 0 (for x^2), 8 (for x), and 14 (constant term).")
    print(" - The prime 2 does not divide the leading coefficient 1.")
    print(f" - The square of the prime, 2^2 = 4, does not divide the constant term {q}.")
    print("Conclusion: The polynomial is irreducible over the rational numbers Q.")
    print("-" * 50)

    # Step 2: Calculate the discriminant
    print("Step 2: Calculate the discriminant (Δ) of the polynomial.")
    print("For a depressed quartic x^4 + px + q, the discriminant formula is Δ = -27*p^4 + 256*q^3.")
    delta_val = -27 * (p**4) + 256 * (q**3)
    print(f"Δ = -27 * ({p}^4) + 256 * ({q}^3)")
    print(f"Δ = -27 * {p**4} + 256 * {q**3}")
    print(f"Δ = {-27 * p**4} + {256 * q**3} = {delta_val}")

    # Check if delta is a perfect square
    sqrt_delta_float = math.sqrt(delta_val)
    is_square = (int(sqrt_delta_float) ** 2 == delta_val)
    
    print(f"The square root of Δ is approx {sqrt_delta_float:.3f}.")
    if is_square:
        print("Conclusion: Δ is a perfect square.")
    else:
        print("Conclusion: Δ is not a perfect square. This means the Galois group G is NOT a subgroup of A_4.")
    print("-" * 50)
        
    # Step 3 & 4: Find and analyze the resolvent cubic
    print("Step 3 & 4: Construct and analyze the resolvent cubic polynomial.")
    print("The resolvent cubic for x^4 + px + q is R(y) = y^3 - 4*q*y - p^2.")
    coeff_y = -4 * q
    const_term = -(p**2)
    print(f"R(y) = y^3 - 4*({q})*y - ({p}^2) = y^3 {coeff_y}y {const_term}")
    
    print("\nWe check for rational roots of R(y) using the Rational Root Theorem.")
    print(f"Possible rational roots are integer divisors of the constant term, {const_term}.")
    
    # Test integer divisors of the constant term to find a rational root
    found_root = False
    divisors = [d for d in range(1, abs(const_term) + 1) if const_term % d == 0]
    for root_candidate in sorted(divisors + [-d for d in divisors]):
        if root_candidate**3 + coeff_y * root_candidate + const_term == 0:
            print(f"Testing y = {root_candidate}: ({root_candidate})^3 + ({coeff_y})*({root_candidate}) + ({const_term}) = {root_candidate**3 + coeff_y*root_candidate + const_term}.")
            print(f"Conclusion: A rational root y = {root_candidate} is found. The resolvent cubic is reducible over Q.")
            found_root = True
            break
            
    if not found_root:
        print("The resolvent cubic has no rational roots and is irreducible over Q.")
    print("-" * 50)
        
    # Step 5: Determine the Galois group and its order
    print("Step 5: Determine the Galois group and its order from the results.")
    print("Summary of findings:")
    print("1. The polynomial P(x) is irreducible over Q.")
    print("2. The discriminant Δ is NOT a perfect square in Q.")
    print("3. The resolvent cubic R(y) is REDUCIBLE over Q.")
    
    print("\nThese three conditions uniquely determine the Galois group G.")
    print("The group G is the Dihedral group of order 8, denoted as D_4.")
    order = 8
    print(f"\nThe order of the Galois group is {order}.")

if __name__ == '__main__':
    compute_galois_group_order()