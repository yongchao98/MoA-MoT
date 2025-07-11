import math

def combinations(n, k):
    """
    Calculates the number of combinations (n choose k).
    Returns 0 if k < 0 or k > n, as is standard.
    """
    if k < 0 or k > n:
        return 0
    # For integer arithmetic, this is more stable than using factorials directly
    # with potential floating point issues, although Python handles large integers.
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def solve():
    """
    Calculates the number of possible score sequences based on generating functions.
    The formula is: Sum_{i,j} C(8,i)C(6,j)(-1)^{i+j} * C(48-7i-9j+13, 13)
    """
    team_a_size = 8
    team_b_size = 6
    
    # max_a_wins = team_b_size
    # max_b_wins = team_a_size
    
    total_games = team_a_size * team_b_size
    
    # Power for the (1-x)^-N term from the denominator
    # N = team_a_size + team_b_size = 8 + 6 = 14
    N = team_a_size + team_b_size
    
    total_sequences = 0
    
    # This equation represents the coefficient of x^48 in the expansion of
    # (1-x^7)^8 * (1-x^9)^6 * (1-x)^(-14).
    # i is the power of x^7, j is the power of x^9.
    
    print("Calculating the total number of sequences.")
    print("The formula is Sum[ C(8,i) * C(6,j) * (-1)^(i+j) * C(k+13, 13) ] where k = 48-7i-9j >= 0.")
    print("---------------------------------------------------------------")
    
    for i in range(team_a_size + 1): # i from 0 to 8
        for j in range(team_b_size + 1): # j from 0 to 6
            
            # The exponent on x from the (1-x^k) terms
            exponent_from_poly = 7 * i + 9 * j
            
            # The required exponent from the (1-x)^-N part
            k = total_games - exponent_from_poly
            
            if k < 0:
                continue

            # Calculate term from (1-x^7)^8
            term1 = combinations(team_a_size, i)
            
            # Calculate term from (1-x^9)^6
            term2 = combinations(team_b_size, j)
            
            # sign from (-1)^i and (-1)^j
            sign = (-1)**(i + j)
            
            # Calculate term from (1-x)^-14, which is C(k+N-1, N-1) = C(k+13, 13)
            term3 = combinations(k + N - 1, N - 1)
            
            term_value = term1 * term2 * sign * term3
            
            if term_value != 0:
                print(f"i={i}, j={j}: C(8,{i})*C(6,{j})*({sign})*C({k+13},{13}) = {term_value}")

            total_sequences += term_value

    print("---------------------------------------------------------------")
    print(f"The total number of different sequences is:")
    print(total_sequences)

solve()
<<<3828353322474496>>>