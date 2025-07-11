import math

def calculate_tori_number_formula():
    """
    This script calculates and presents the formula for the exact number of
    Fq-rational maximal tori in a reductive group of type E_8.
    """

    # 1. Define the parameters for the Lie group of type E_8.
    # The rank of E8 is r.
    rank = 8
    # The degrees of the fundamental invariants of the Weyl group W(E_8).
    degrees = [2, 8, 12, 14, 18, 20, 24, 30]
    # The number of positive roots of E_8.
    # It is calculated as (dim(G) - rank) / 2 = (248 - 8) / 2.
    num_positive_roots = 120

    # 2. Calculate the order of the Weyl group W(E_8).
    # The order of a Weyl group is the product of the degrees of its fundamental invariants.
    weyl_order = math.prod(degrees)

    # 3. Construct the final formula string with all numbers explicitly included.
    # The numerator is |G(F_q)| = q^N * Product[i=1 to r](q^d_i - 1).
    numerator_terms = ' * '.join(f"(q^{d} - 1)" for d in degrees)
    
    # The denominator is |N_G(T_0)(F_q)| = |W(E_8)| * (q - 1)^r.
    
    # The full formula is N(q) = |G(F_q)| / |N_G(T_0)(F_q)|.
    final_formula = f"N(q) = (q^{num_positive_roots} * {numerator_terms}) / ({weyl_order} * (q - 1)^{rank})"
    
    print("The formula for the exact number of rational maximal tori of a group of type E_8 over F_q is:")
    print(final_formula)

if __name__ == '__main__':
    calculate_tori_number_formula()