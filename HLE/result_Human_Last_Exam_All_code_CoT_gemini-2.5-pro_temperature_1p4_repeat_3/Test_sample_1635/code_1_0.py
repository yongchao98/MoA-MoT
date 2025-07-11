def solve_cardinality_S():
    """
    Solves the problem based on Sharkovsky's Theorem and number theory.
    """
    # Step 1: From Sharkovsky's theorem, we deduce the set of non-existent least periods.
    # The existence of a point of period 13 and non-existence of a point of period 11
    # implies that the only least periods that DO NOT exist are 3, 5, 7, 9, and 11.
    S_least = {3, 5, 7, 9, 11}

    # Step 2: Define the set S based on the problem statement.
    # k is in S if and only if for any point x, f^k(x)=x implies f(x)=x.
    # This means that no divisor of k (other than 1) can be an existing least period.
    # So, all divisors of k (other than 1) must be in S_least.
    
    # We will search for such numbers k.
    # A number k > 11 cannot be in S, because k itself is a divisor of k,
    # and k would not be in S_least. So we only need to search up to k=11.
    # We search a bit beyond to be fully confident.
    
    S = []
    search_limit = 100 
    
    for k in range(1, search_limit):
        # Find all divisors of k greater than 1
        divisors_of_k = set()
        for i in range(2, k + 1):
            if k % i == 0:
                divisors_of_k.add(i)
        
        # Check if all these divisors are in S_least
        if divisors_of_k.issubset(S_least):
            S.append(k)
            
    print("The set S = {k : there is no point of order k} is:")
    # Using print to output each number in the set as requested
    print(S)
    
    print("\nThe cardinality of S is:")
    print(len(S))

solve_cardinality_S()
<<<6>>>