import math

def calculate_dimension(n):
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1(2) on P^n_C.
    """
    if n < 0:
        print("n must be a non-negative integer.")
        return

    # h^0(P^n, O(1)^{n+1}) = (n+1) * h^0(P^n, O(1)) = (n+1) * (n+1)
    h0_O1_sum = (n + 1)**2

    # h^0(P^n, O(2)) = C(n+2, 2)
    h0_O2 = math.comb(n + 2, 2)

    # The dimension is the difference due to the short exact sequence of global sections.
    dimension = h0_O1_sum - h0_O2

    print(f"To find the dimension for n = {n}, we use the formula:")
    print(f"dim = h0(P^{n}, O(1)^{{n+1}}) - h0(P^{n}, O(2))")
    print(f"dim = (n+1)^2 - C(n+2, 2)")
    print(f"\nCalculating the terms for n = {n}:")
    print(f"h0(P^{n}, O(1)^{{n+1}}) = ({n}+1)^2 = {h0_O1_sum}")
    print(f"h0(P^{n}, O(2)) = C({n}+2, 2) = {h0_O2}")
    
    print("\nThe final equation with these numbers is:")
    print(f"{h0_O1_sum} - {h0_O2} = {dimension}")

    print(f"\nThe general formula for the dimension is n(n+1)/2, which is C(n+1, 2).")
    print(f"For n={n}, this is C({n+1}, 2) = {math.comb(n+1, 2)}, which matches our result.")


# We can choose any integer n >= 0. Let's use n=4 as an example.
n_value = 4
calculate_dimension(n_value)
