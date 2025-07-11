import math

def calculate_dimension(n):
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1(2) on P^n_C.

    The dimension is given by the formula: h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2)).
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # h^0(P^n, O(k)) = C(n+k, k). We can use polynomial forms for simplicity.
    # h^0(P^n, O(1)) = n + 1
    h0_O1 = n + 1
    
    # h^0(P^n, O(2)) = (n+2)*(n+1)/2
    h0_O2 = (n + 2) * (n + 1) // 2
    
    # h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1))
    term1 = (n + 1) * h0_O1
    
    # The dimension is the difference
    dimension = term1 - h0_O2
    
    print(f"For n = {n} (the complex projective space P^{n}):")
    print("The dimension of H^0(P^n, Omega^1(2)) is calculated from the exact sequence of global sections.")
    print("dim = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))\n")
    
    print("First term calculation:")
    print(f"h^0(P^n, O(1)) = n+1 = {h0_O1}")
    print(f"h^0(P^n, O(1)^(n+1)) = (n+1) * {h0_O1} = {term1}\n")
    
    print("Second term calculation:")
    print(f"h^0(P^n, O(2)) = (n+2)(n+1)/2 = {h0_O2}\n")
    
    print("Final equation with calculated values:")
    print(f"dim = {term1} - {h0_O2} = {dimension}")
    
    # The formula simplifies to n(n+1)/2
    simplified_dimension = n * (n + 1) // 2
    assert dimension == simplified_dimension

# The problem is posed for a general integer n.
# As a concrete example, we calculate the dimension for the complex projective plane, n=2.
n_value = 2
calculate_dimension(n_value)
