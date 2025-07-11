import math
from fractions import Fraction

def solve_mass_formula():
    """
    This function calculates the total mass M(A_5, rho, 2)
    based on established mass formulas in number theory.
    """
    
    # Define parameters from the problem
    G_name = "A_5"
    n = 5
    p = 2
    
    # The order of the group G = A_5
    order_G = math.factorial(n) // 2
    
    # The order of the symmetric group S_5
    order_Sn = math.factorial(n)

    # Step 1: Use Serre's mass formula for the S_n permutation representation.
    # The sum of |Aut(E)|*|disc(E)|_p over all quintic etale algebras E is |S_5|.
    S_total = order_Sn
    
    # Step 2: An etale algebra is an A_n-algebra if and only if its discriminant is a square.
    # We need to know the number of square classes in Q_2.
    # For a p-adic field Q_p, the group Q_p*/(Q_p*)^2 has order 8 for p=2.
    num_square_classes_in_Q2 = 8
    
    # Step 3: A theorem from local class field theory implies that the sum S_total
    # is equidistributed among the square classes. So, to get the sum for A_5 algebras,
    # we take the total sum and divide by the number of square classes.
    S_A5 = S_total / num_square_classes_in_Q2
    
    # Step 4: The mass M(A_5, rho, 2) is this sum for A_5 algebras divided by the order of A_5.
    mass = Fraction(int(S_A5), order_G)

    print("The total mass M(A_5, rho, 2) is computed by the following steps:")
    print("1. The sum over all quintic etale algebras, Sum_total = |S_5| = 120.")
    print(f"2. The number of square classes in Q_2 is {num_square_classes_in_Q2}.")
    print("3. The sum for A_5-algebras (square discriminant class) is Sum_A5 = Sum_total / number of square classes.")
    print("4. The order of A_5 is |A_5| = 60.")

    print("\nThe final calculation is:")
    # The final equation outputting each number.
    print(f"M(A_5, rho, 2) = {int(S_A5)} / {order_G} = {mass.numerator}/{mass.denominator}")

solve_mass_formula()