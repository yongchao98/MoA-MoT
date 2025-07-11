import sympy
from sympy import SymmetricGroup, Rational
from sympy.combinatorics import Permutation

def solve_hurwitz_number():
    """
    Calculates the Hurwitz number for a degree 7 cover of P^1 with four specified
    ramification profiles using the character theory of the symmetric group.
    """
    d = 7
    G = SymmetricGroup(d)

    # Define permutations with the required cycle structures to represent the conjugacy classes
    # Ramification profiles: (4,3), (5,2), (2,1^5), (2,1^5)
    perm1 = Permutation(cycles=[[0, 1, 2, 3], [4, 5, 6]])  # Cycle type (4,3)
    perm2 = Permutation(cycles=[[0, 1, 2, 3, 4], [5, 6]])   # Cycle type (5,2)
    perm3 = Permutation(cycles=[[0, 1]])                  # Cycle type (2)

    # Get the sizes of the corresponding conjugacy classes
    size_C1 = G.get_conjugacy_class(perm1).size
    size_C2 = G.get_conjugacy_class(perm2).size
    size_C3 = G.get_conjugacy_class(perm3).size
    # The fourth class is the same as the third
    size_C4 = size_C3
    
    # Get the order of the group
    order_G = G.order()

    # Print the formula and the calculated components
    print("The formula for the Hurwitz number N is:")
    print("N = (|C(4,3)| * |C(5,2)| * |C(2)| * |C(2)| / |S_7|) * Sum")
    print("where Sum = sum over all irreducible characters chi of (chi(C(4,3)) * chi(C(5,2)) * chi(C(2))**2) / chi(1)**2\n")

    print("Calculating the values for the formula:")
    print(f"Order of the group S_7: |S_7| = {order_G}")
    print(f"Size of conjugacy class for partition (4,3): |C(4,3)| = {size_C1}")
    print(f"Size of conjugacy class for partition (5,2): |C(5,2)| = {size_C2}")
    print(f"Size of conjugacy class for partition (2,1,1,1,1,1): |C(2)| = {size_C3}\n")

    prefactor_val = Rational(size_C1 * size_C2 * size_C3 * size_C4, order_G)
    print("The prefactor in the formula is:")
    print(f"(|C(4,3)| * |C(5,2)| * |C(2)|**2) / |S_7| = ({size_C1} * {size_C2} * {size_C3}**2) / {order_G} = {prefactor_val}\n")

    # Calculate the sum over all irreducible characters
    char_table = G.character_table()
    total_sum = Rational(0)

    for chi in char_table:
        # chi(perm) gives the character value for the class of perm
        # chi(G.identity) gives chi(1), the character dimension
        chi_1 = chi(G.identity)
        chi_C1 = chi(perm1)
        chi_C2 = chi(perm2)
        chi_C3 = chi(perm3)
        
        numerator = chi_C1 * chi_C2 * chi_C3**2
        denominator = chi_1**2
        
        term = Rational(numerator, denominator)
        total_sum += term

    print("Calculating the sum over characters:")
    print(f"Sum = {total_sum.p}/{total_sum.q} = {float(total_sum)}\n")

    # Final calculation for the Hurwitz number
    hurwitz_number = prefactor_val * total_sum
    
    print("The final Hurwitz number is:")
    print(f"N = ({prefactor_val.p}/{prefactor_val.q}) * ({total_sum.p}/{total_sum.q})")
    print(f"N = {hurwitz_number}")

solve_hurwitz_number()