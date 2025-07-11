import math

def combinations(n, k):
    """
    Calculates the number of combinations 'n choose k', C(n,k).
    This implementation avoids large intermediate numbers and floating point errors.
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1(2) on P^n_C, and explains the steps.
    """
    try:
        n_str = input("Please enter the dimension n of the complex projective space (a non-negative integer): ")
        n = int(n_str)
        if n < 0:
            print("Error: The dimension n must be a non-negative integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    print(f"\nWe want to compute h^0(P^{n}, Omega^1(2)), the dimension of the space of global sections.")
    print("This can be found using the Euler sequence and sheaf cohomology.")
    print("The key relation is derived from a short exact sequence of vector spaces, which gives the formula:")
    print(f"h^0(P^{n}, Omega^1(2)) = h^0(P^{n}, O(1)^(n+1)) - h^0(P^{n}, O(2))\n")
    
    print(f"For n = {n}:")
    
    # h^0(P^n, O(1)^(n+1))
    h0_O1 = combinations(n + 1, 1)
    h0_O1_n1 = (n + 1) * h0_O1

    # h^0(P^n, O(2))
    h0_O2 = combinations(n + 2, 2)
    
    # Result
    result = h0_O1_n1 - h0_O2
    
    print("1. Calculate the dimension of H^0(P^n, O(1)^(n+1)):")
    print(f"   h^0(P^{n}, O(1)) = C(n+1, 1) = C({n+1}, 1) = {h0_O1}")
    print(f"   h^0(P^{n}, O(1)^(n+1)) = (n+1) * h^0(P^{n}, O(1)) = {n+1} * {h0_O1} = {h0_O1_n1}")
    
    print("\n2. Calculate the dimension of H^0(P^n, O(2)):")
    print(f"   h^0(P^{n}, O(2)) = C(n+2, 2) = C({n+2}, 2) = {h0_O2}")
    
    print("\n3. Subtract the second value from the first to get the final dimension:")
    print(f"   h^0(P^{n}, Omega^1(2)) = {h0_O1_n1} - {h0_O2} = {result}")

    print(f"\nIn general, for a given n, the dimension is given by the formula n(n+1)/2, which is C(n+1, 2).")

if __name__ == '__main__':
    calculate_dimension()