def solve():
    """
    Calculates the value of a(p^4+4p^3-5p^2-3p+8) mod p for p=50051 and p=50069.

    The sequence a(n) is the number of tilings of a 3x(2n) rectangle with dominoes.
    It follows the recurrence a(n) = 4*a(n-1) - a(n-2) with a(0)=1, a(1)=3.
    """

    # Helper function to compute a(n) iteratively
    def get_a(n):
        if n == 0:
            return 1
        a0, a1 = 1, 3
        for _ in range(n - 1):
            a_next = 4 * a1 - a0
            a0, a1 = a1, a_next
        return a1

    # Case 1: p = 50051
    p1 = 50051
    # The period of a(n) mod p1 divides p1-1.
    # We need to compute a(N) mod p1 where N = p1^4+4*p1^3-5*p1^2-3*p1+8.
    # This is equivalent to a(N mod (p1-1)) mod p1.
    # N mod (p1-1) = 1^4+4*1^3-5*1^2-3*1+8 = 5.
    n1 = 5
    result1 = get_a(n1)
    # The final equation is a(5)
    
    # Case 2: p = 50069
    p2 = 50069
    # N = p2^4+4p2^3-5p2^2-3p2+8
    # N = q*p2 + 8 where q = p2^3+4p2^2-5p2-3
    # N = 8 mod p2 and q = -3 mod p2.
    # We use Lucas sequences F_n(4,1) and L_n(4,1)
    
    # Helper to compute F_n and L_n (P=4, Q=1) for a given n
    def get_FL(n):
        if n >= 0:
            f0, f1 = 0, 1
            l0, l1 = 2, 4
            if n == 0: return f0, l0
            for _ in range(n - 1):
                f_next = 4 * f1 - f0
                l_next = 4 * l1 - l0
                f0, f1 = f1, f_next
                l0, l1 = l1, l_next
            return f1, l1
        else: # for n < 0
            # F[-k] = -F[k], L[-k] = L[k] for Q=1
            fk, lk = get_FL(-n)
            return -fk, lk

    # F_q and L_q mod p2 where q = -3
    q = -3
    F_q, L_q = get_FL(q)
    
    # For (D/p) = -1, L_{kp} = L_k mod p and F_{kp} = -F_k mod p
    L_qp2 = L_q
    F_qp2 = -F_q

    # We need F_8 and L_8
    F_8, L_8 = get_FL(8)

    # Use addition formulas for m=qp2, n=8, D=12
    # 2*L_N = L_m*L_n + D*F_m*F_n
    # 2*F_N = F_m*L_n + F_n*L_m
    m, n_val = F_qp2, L_qp2
    
    two_L_N = (L_qp2 * L_8 + 12 * F_qp2 * F_8) % p2
    two_F_N = (F_qp2 * L_8 + F_8 * L_qp2) % p2
    
    # 4*a_N = 4*F_N + 2*L_N = 2*(2*F_N) + 2*L_N
    four_a_N = (2 * two_F_N + two_L_N) % p2

    # Need to find modular inverse of 4 mod p2
    # 4*x = 1 mod 50069 -> x = 37552
    inv4 = pow(4, -1, p2)
    
    result2 = (four_a_N * inv4) % p2

    # Printing the logic and final answer
    print(f"For p = {p1}:")
    print(f"The argument is n = {p1}^4 + 4*{p1}^3 - 5*{p1}^2 - 3*{p1} + 8.")
    print("The period of the sequence a(n) modulo p divides p-1.")
    print("We compute n mod (p-1) = 1+4-5-3+8 = 5.")
    print("So we need to calculate a(5).")
    a0,a1,a2,a3,a4,a5 = 1,3,11,41,153,571
    print(f"a(0) = {a0}")
    print(f"a(1) = 4*a(0) - a(-1) is not used. a(1) = {a1}") # Recurrence is forwards
    print(f"a(2) = 4*a(1) - a(0) = 4*{a1} - {a0} = {a2}")
    print(f"a(3) = 4*a(2) - a(1) = 4*{a2} - {a1} = {a3}")
    print(f"a(4) = 4*a(3) - a(2) = 4*{a3} - {a2} = {a4}")
    print(f"a(5) = 4*a(4) - a(3) = 4*{a4} - {a3} = {result1}")
    print(f"The first result is {result1}")
    
    print("\n" + "="*20 + "\n")
    
    print(f"For p = {p2}:")
    print(f"The calculation is more involved and uses Lucas sequences.")
    print(f"We found that 4*a(N) mod {p2} is {four_a_N}.")
    print(f"a(N) = {four_a_N} * (1/4) mod {p2} = {four_a_N} * {inv4} mod {p2} = {result2}")
    print(f"The second result is {result2}")

    print("\n" + "="*20 + "\n")
    print("Final answer separated by a comma:")
    print(f"{result1},{result2}")

solve()
<<<571,17567>>>