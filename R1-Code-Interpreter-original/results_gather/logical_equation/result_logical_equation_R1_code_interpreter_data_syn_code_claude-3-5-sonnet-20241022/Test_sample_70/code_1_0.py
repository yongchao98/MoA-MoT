def find_solution():
    # Available numbers
    numbers = set([1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50])
    result = [0] * 13  # A to M
    used = set()

    # From equation 6: L + M = 51
    # And equation 1: M = 2.4L
    # We can solve: L + 2.4L = 51
    # Therefore: 3.4L = 51
    # L = 15
    L = 15
    M = 36  # M = 2.4 * 15 = 36
    result[11] = L  # L
    result[12] = M  # M
    used.add(L)
    used.add(M)

    # From equation 2: F = 3.0L
    F = int(3.0 * L)  # F = 45
    result[5] = F
    used.add(F)

    # From equation 10: H = 1.5K
    # And equation 4: K - H = -8
    # Substituting: K - 1.5K = -8
    # -0.5K = -8
    # K = 16
    K = 16
    H = int(1.5 * K)  # H = 24
    result[7] = H
    result[10] = K
    used.add(K)
    used.add(H)

    # From equation 5: H > I
    # And equation 3: E + I = 8
    # Since H = 24, I must be < 24
    # Given I + E = 8, and both must be from our number set
    # I = 5, E = 3 is the only possibility
    I = 5
    E = 3
    result[4] = E
    result[8] = I
    used.add(E)
    used.add(I)

    # From equation 8: G = 2.0I
    G = int(2.0 * I)  # G = 10
    result[6] = G
    used.add(G)

    # From equation 7: C + J = 9
    # Looking at remaining numbers, only 2 + 7 = 9
    C = 2
    J = 7
    result[2] = C
    result[9] = J
    used.add(C)
    used.add(J)

    # From equation 9: A - D = 49
    # Looking at remaining numbers: 50 - 1 = 49
    A = 50
    D = 1
    result[0] = A
    result[3] = D
    used.add(A)
    used.add(D)

    # B is the only remaining number
    B = (numbers - used).pop()
    result[1] = B

    print(result)

find_solution()