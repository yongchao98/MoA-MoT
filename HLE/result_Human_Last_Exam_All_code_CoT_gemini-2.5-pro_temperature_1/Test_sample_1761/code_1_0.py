import math

def combinations(n, k):
    """
    Calculates n C k, the number of combinations.
    """
    if k < 0 or k > n:
        return 0
    # Use math.comb for Python 3.8+ for efficiency and safety
    if hasattr(math, 'comb'):
        return math.comb(n, k)
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_dimension():
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1(2) on P^n, the complex projective space.
    """
    try:
        n_str = input("Enter the dimension n of the complex projective space P^n (e.g., n=2 for the projective plane): ")
        n = int(n_str)
        if n < 1:
            print("Error: The dimension n must be an integer greater than or equal to 1.")
            return
    except ValueError:
        print("Error: Please enter a valid integer for n.")
        return

    # h^0(P^n, O(1)) = C(n+1, 1) = n+1
    h0_O1 = combinations(n + 1, 1)
    
    # h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1))
    term1 = (n + 1) * h0_O1
    
    # h^0(P^n, O(2)) = C(n+2, 2)
    term2 = combinations(n + 2, 2)
    
    # The result is the difference
    result = term1 - term2
    
    print(f"\nFor n = {n}, the dimension is calculated using the formula:")
    print(f"h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))\n")
    
    print(f"First term: h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1))")
    print(f"                                = ({n}+1) * C({n}+1, 1)")
    print(f"                                = {n+1} * {h0_O1} = {term1}")
    
    print(f"\nSecond term: h^0(P^n, O(2)) = C(n+2, 2) = C({n+2}, 2) = {term2}")
    
    print("\nFinal equation:")
    print(f"h^0(P^n, Omega^1(2)) = {term1} - {term2} = {result}\n")
    
    # Also print the general formula result for verification
    formula_result = n * (n + 1) // 2
    print(f"This matches the general simplified formula for the dimension: n(n+1)/2 = {n}({n}+1)/2 = {formula_result}")

# Execute the function
solve_dimension()