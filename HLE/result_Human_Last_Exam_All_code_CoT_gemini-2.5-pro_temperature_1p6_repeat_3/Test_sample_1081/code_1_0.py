import math

def solve_e8_tori():
    """
    This script calculates the components of the formula for the number of
    Fq-rational maximal tori of a reductive group of type E8.
    """

    # For a group of type E8, the rank is 8.
    rank = 8
    
    # The dimension of the E8 group is 248.
    # The number of roots is dim(G) - rank = 248 - 8 = 240.
    # The number of positive roots N is half the number of roots.
    num_positive_roots = 120
    
    # The degrees of the fundamental invariants for E8.
    degrees = [2, 8, 12, 14, 18, 20, 24, 30]
    
    # The order of the Weyl group for E8 is the product of the degrees.
    weyl_order = math.prod(degrees)
    
    print("The problem is to find the number of rational maximal tori of a group G of type E8 over a finite field Fq.")
    print("The group G is necessarily the split form of E8.")
    print("\nThe key components of the formula are:")
    print(f"1. Rank of the group (r): {rank}")
    print(f"2. Number of positive roots (N): {num_positive_roots}")
    print(f"3. Degrees of fundamental invariants (d_i): {degrees}")
    print(f"4. Order of the Weyl group W(E8): {weyl_order}")

    print("\nThe number of rational maximal tori is given by the formula |G(Fq)| / |N(T0)(Fq)|.")
    
    # Constructing the formula string
    # Numerator: |G(Fq)| = q^N * product(q^d_i - 1)
    numerator_str_parts = [f"(q^{d} - 1)" for d in degrees]
    numerator_str = f"q^{num_positive_roots} * " + " * ".join(numerator_str_parts)
    
    # Denominator: |N(T0)(Fq)| = (q-1)^r * |W(E8)|
    denominator_str = f"(q - 1)^{rank} * {weyl_order}"
    
    print("\nFinal formula for the number of tori:")
    print(f"N = ({numerator_str}) / ({denominator_str})")

solve_e8_tori()