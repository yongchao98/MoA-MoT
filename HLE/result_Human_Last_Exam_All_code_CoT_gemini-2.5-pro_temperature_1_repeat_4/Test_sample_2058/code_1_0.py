from fractions import Fraction

def solve_mass_problem():
    """
    Calculates the total mass M(A_5, rho, 2) based on mass formulas.
    """
    
    # Order of the group G = A_5
    order_A5 = 60
    
    # Order of the symmetric group S_5
    order_S5 = 120
    
    print(f"The order of A_5 is |A_5| = {order_A5}.")
    print(f"The order of S_5 is |S_5| = {order_S5}.")
    print("-" * 20)

    # The sum over all degree 5 etale algebras is |S_5|.
    # Let S(S_5) be this sum.
    S_S5 = order_S5
    print(f"The total sum over homomorphisms from Gamma_2 to S_5 is S(S_5) = {S_S5}.")

    # Number of square classes in Q_2^* / (Q_2^*)^2 is 8.
    num_square_classes_Q2 = 8
    print(f"The number of square classes in Q_2 is {num_square_classes_Q2}.")
    
    # By the principle of equidistribution of discriminants, the sum for maps
    # into A_5 (where the discriminant is a square) is the total sum divided
    # by the number of square classes.
    # Let S(A_5) be this sum.
    S_A5 = S_S5 / num_square_classes_Q2
    
    print(f"The sum S(A_5) for homomorphisms into A_5 is S(S_5) / {num_square_classes_Q2}.")
    print(f"S(A_5) = {S_S5} / {num_square_classes_Q2} = {int(S_A5)}.")
    print("-" * 20)

    # The total mass is M = S(A_5) / |A_5|.
    mass = Fraction(int(S_A5), order_A5)
    
    print("The total mass M(A_5, rho, 2) is S(A_5) / |A_5|.")
    print(f"M(A_5, rho, 2) = {int(S_A5)} / {order_A5} = {mass.numerator}/{mass.denominator}")

solve_mass_problem()