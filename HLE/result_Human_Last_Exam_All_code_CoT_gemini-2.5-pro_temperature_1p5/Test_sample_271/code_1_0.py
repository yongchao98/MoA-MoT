import sys

def solve():
    """
    Solves the problem of finding all possible degrees of intermediate normal
    field extensions for the splitting field of P(x) over Q.
    """
    print("Step 1: Analyze the polynomial P(x) = x^7 - 2x^5 - 9x^3 + 3x^2 + 18x - 6.")
    print("By inspection, we can factor the polynomial. Let's verify the factorization:")
    print("P(x) = (x^5 - 9x + 3) * (x^2 - 2)")
    print("Expanding this gives: x^7 - 2x^5 - 9x^3 + 18x + 3x^2 - 6, which is P(x).")
    print("So, the roots of P(x) are the roots of f(x) = x^5 - 9x + 3 and g(x) = x^2 - 2.")
    print("-" * 20)

    print("Step 2: Determine the Galois group G = Gal(K/Q).")
    print("The splitting field K is the compositum of the splitting fields of f(x) and g(x).")
    
    print("\nAnalysis of g(x) = x^2 - 2:")
    print("The roots are +/- sqrt(2). The splitting field is Kg = Q(sqrt(2)).")
    print("The Galois group Gg = Gal(Kg/Q) is isomorphic to S_2 (or Z/2Z), with order 2.")
    
    print("\nAnalysis of f(x) = x^5 - 9x + 3:")
    print("1. Irreducibility: By Eisenstein's criterion with prime p=3, f(x) is irreducible over Q.")
    print("2. Number of real roots: f'(x) = 5x^4 - 9 has two real roots, so f(x) has at most 3 real roots. Evaluating f(x) at its critical points shows it has exactly 3 real roots and 2 complex roots.")
    print("An irreducible quintic polynomial over Q with exactly 3 real roots has a Galois group isomorphic to S_5.")
    print("So, the Galois group Gf = Gal(Kf/Q) is S_5, which has order 5! = 120.")

    print("\nCombining the Galois groups:")
    print("The Galois group G = Gal(K/Q) is a subgroup of Gf x Gg = S_5 x S_2.")
    print("If the intersection of the splitting fields Kf and Kg is just Q, then G is the direct product S_5 x S_2.")
    print("The intersection Kf intersect Kg is non-trivial if and only if Kg is a subfield of Kf.")
    print("Kg = Q(sqrt(2)) is the unique quadratic subfield of Kf if and only if the discriminant of f(x) is of the form c^2 * 2 for some c in Q.")
    print("Let's calculate the discriminant of f(x) = x^5 + ax + b (with a=-9, b=3).")
    n, a, b = 5, -9, 3
    # Formula for discriminant of x^n+ax+b is (-1)^(n(n-1)/2) * (n^n * b^(n-1) + (-1)^(n-1) * (n-1)^(n-1) * a^n)
    disc = ((-1)**(n*(n-1)//2)) * (n**n * b**(n-1) + ((-1)**(n-1)) * (n-1)**(n-1) * a**n)
    print(f"The discriminant of f(x) is {disc}.")
    print("This is a negative integer, so sqrt(disc) is imaginary. Thus, the quadratic subfield of Kf is Q(sqrt(disc)), which is an imaginary quadratic field.")
    print("Since Q(sqrt(2)) is a real field, it cannot be the quadratic subfield of Kf.")
    print("Therefore, Kf intersect Kg = Q, and the Galois group G = Gal(K/Q) is isomorphic to S_5 x S_2.")
    order_G = 120 * 2
    print(f"The order of G is 120 * 2 = {order_G}.")
    print("-" * 20)

    print("Step 3: Find normal subgroups of G and their indices.")
    print("By the Fundamental Theorem of Galois Theory, normal extensions L/Q with Q < L < K correspond to proper non-trivial normal subgroups H of G.")
    print("The degree of the extension [L:Q] is the index of the subgroup H in G, i.e., [G:H] = |G|/|H|.")
    print("The normal subgroups of S_5 are {e}, A_5, S_5 (orders 1, 60, 120).")
    print("The normal subgroups of S_2 are {e}, S_2 (orders 1, 2).")
    
    degrees = set()
    
    print("\nCase 1: Normal subgroups that are direct products of normal subgroups of S_5 and S_2.")
    # H = N1 x N2, where N1 is normal in S5 and N2 is normal in S2.
    # We exclude H = {e}x{e} (gives K) and H = S5xS2 (gives Q).
    # H_1 = {e} x S_2: |H|=2, index=240/2=120. Corresponds to L=Kf.
    degrees.add(order_G // 2)
    print("H = {e} x S_2 -> |H|=2. Degree = 240/2 = 120.")
    # H_2 = A_5 x {e}: |H|=60, index=240/60=4.
    degrees.add(order_G // 60)
    print("H = A_5 x {e} -> |H|=60. Degree = 240/60 = 4.")
    # H_3 = A_5 x S_2: |H|=120, index=240/120=2.
    degrees.add(order_G // 120)
    print("H = A_5 x S_2 -> |H|=120. Degree = 240/120 = 2.")
    # H_4 = S_5 x {e}: |H|=120, index=240/120=2. Corresponds to L=Kg.
    degrees.add(order_G // 120)
    print("H = S_5 x {e} -> |H|=120. Degree = 240/120 = 2.")

    print("\nCase 2: 'Diagonal' normal subgroups.")
    print("There can be other normal subgroups when there's an isomorphism between quotient groups of factors. Here, S_5/A_5 is isomorphic to S_2/{e} (both are S_2).")
    print("This gives rise to a normal subgroup H = {(g,s) in S_5 x S_2 | sgn(g)=s}. The order of this subgroup is 120.")
    # H_5 has |A_5| elements where s=e and |S_5 \ A_5| elements where s is not e. So |H|=60+60=120.
    degrees.add(order_G // 120)
    print("H_diag -> |H|=120. Degree = 240/120 = 2.")
    print("-" * 20)

    print("Step 4: List the unique possible degrees.")
    # Sort the degrees for a clean presentation.
    sorted_degrees = sorted(list(degrees))
    print("The set of all possible degrees for such extensions L is the set of these indices.")
    print("The possible degrees are:", ", ".join(map(str, sorted_degrees)))

if __name__ == '__main__':
    solve()
