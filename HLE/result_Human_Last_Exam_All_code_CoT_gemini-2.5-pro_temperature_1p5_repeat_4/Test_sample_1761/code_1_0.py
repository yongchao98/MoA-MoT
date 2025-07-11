import math

def solve_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_{P^n} tensored with O(2) on complex projective space P^n.
    """
    # The dimension 'n' of the complex projective space P^n.
    # The problem is posed for a general n. We will use n=2 as an example.
    # You can change this value to any integer n >= 1.
    n = 2

    # The dimension is given by the formula: h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2)).
    
    # Calculate h^0(O(1)^(n+1)) which is (n+1)^2.
    h0_O1_oplus_n1 = (n + 1)**2
    
    # Calculate h^0(O(2)) which is C(n+2, 2).
    # math.comb(k, k) calculates the binomial coefficient "n choose k".
    h0_O2 = math.comb(n + 2, 2)
    
    # Calculate the final dimension.
    dimension = h0_O1_oplus_n1 - h0_O2

    print(f"The calculation is for the complex projective space P^n where n = {n}.")
    print("The dimension is calculated using the formula:")
    print("dim = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))")
    print(f"h^0(P^{n}, O(1)^(n+1)) = ({n}+1)^2 = {h0_O1_oplus_n1}")
    print(f"h^0(P^{n}, O(2)) = C({n}+2, 2) = {h0_O2}")
    print("\nThe resulting dimension is:")
    print(f"dim = {h0_O1_oplus_n1} - {h0_O2} = {dimension}")
    
    # The simplified formula is C(n+1, 2)
    simple_dim = math.comb(n+1, 2)
    print(f"\nVerification with the simplified formula C(n+1, 2) = C({n}+1, 2) = {simple_dim}")


solve_dimension()