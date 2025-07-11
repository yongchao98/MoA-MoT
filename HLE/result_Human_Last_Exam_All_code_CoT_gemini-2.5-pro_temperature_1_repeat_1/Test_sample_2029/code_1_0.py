import math

def get_comb():
    """Returns the math.comb function, with a fallback for older Python versions."""
    try:
        # math.comb is available in Python 3.8+ and is efficient.
        from math import comb
        return comb
    except ImportError:
        # Fallback for Python versions older than 3.8.
        def comb_fallback(n, k):
            if k < 0 or k > n:
                return 0
            if k == 0 or k == n:
                return 1
            if k > n // 2:
                k = n - k
            
            # Using integer division
            res = 1
            for i in range(k):
                res = res * (n - i) // (i + 1)
            return res
        return comb_fallback

def solve_match_sequences():
    """
    Calculates the number of possible score sequences in the A vs B match.
    """
    comb = get_comb()

    # Define the problem parameters based on the problem description.
    team_a_members = 8
    team_b_members = 6

    # Max wins for a player is the number of games they play.
    max_wins_a = team_b_members
    max_wins_b = team_a_members

    # The total number of games is the sum of all wins.
    total_games = team_a_members * team_b_members

    # The number of sequences is the coefficient of x^total_games in the expansion of:
    # (1+..+x^max_wins_a)^team_a_members * (1+..+x^max_wins_b)^team_b_members
    # which simplifies to:
    # (1-x^(max_wins_a+1))^team_a_members * (1-x^(max_wins_b+1))^team_b_members * (1-x)^-(team_a_members+team_b_members)

    # Let's define the parameters for the formula
    N1 = team_a_members   # 8
    M1 = max_wins_a + 1    # 7
    N2 = team_b_members   # 6
    M2 = max_wins_b + 1    # 9
    K = total_games        # 48
    EXP = N1 + N2          # 14

    print("The number of sequences is the coefficient of x^48 in the expansion of (1-x^7)^8 * (1-x^9)^6 * (1-x)^-14.")
    print("This is found by summing terms based on the binomial expansions:")
    print("-" * 60)
    
    total_sequences = 0
    
    # Iterate through the expansion of (1-x^M1)^N1
    for i in range(N1 + 1):
        # Iterate through the expansion of (1-x^M2)^N2
        for j in range(N2 + 1):
            
            # The power of x needed from the (1-x)^-EXP expansion
            k = K - M1 * i - M2 * j
            
            if k >= 0:
                # Calculate the coefficient for this combination of i, j, k
                # C(N1, i)*(-1)^i is from (1-x^M1)^N1
                # C(N2, j)*(-1)^j is from (1-x^M2)^N2
                # C(k+EXP-1, EXP-1) is from (1-x)^-EXP
                term_value = (comb(N1, i) * ((-1)**i) *
                              comb(N2, j) * ((-1)**j) *
                              comb(k + EXP - 1, EXP - 1))
                
                # Outputting each number in the final equation's terms
                if term_value != 0:
                    c1_str = f"C({N1},{i})*(-1)^{i}"
                    c2_str = f"C({N2},{j})*(-1)^{j}"
                    c3_str = f"C({k}+{EXP-1},{EXP-1})"
                    print(f"Term for (i={i}, j={j}): [{c1_str}] * [{c2_str}] * [{c3_str}] = {term_value}")

                total_sequences += term_value

    print("-" * 60)
    print("The total number of different sequences is the sum of all terms:")
    print(total_sequences)

solve_match_sequences()