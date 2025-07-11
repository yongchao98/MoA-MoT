import sympy

def solve():
    """
    This function formulates and prints the derived lower bound for m.
    """
    N = sympy.Symbol('N')
    d_prime = sympy.Symbol("d'")

    # Based on the derivation, the lower bound for m is Omega(N * d').
    # We showed that for a fixed set of 'y' inputs, the network must approximate a linear map
    # of rank N*d'. Known results from approximation theory state that the number of hidden
    # neurons 'm' must be at least of the order of the rank of this linear map.
    lower_bound = N * d_prime

    print("The lower bound for the number of neurons m is Omega(N * d').")
    # To express this using the variables:
    # print(f"m >= c * {lower_bound} for some constant c > 0")
    
solve()
