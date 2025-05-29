def verify_solution(solution):
    A, B, C, D, E, F, G, H, I, J, K, L, M = solution
    conditions = [
        A - E == 4,
        E - A == -4,
        I == 4.0 * J,
        D - B == -14,
        H + M == 26,
        abs(L - 3.2 * A) < 0.1,
        abs(L - 1.6 * F) < 0.1,
        F + L == 26,
        H - I == -4,
        abs(D - 1.5 * H) < 0.1
    ]
    return all(conditions)

# Let's solve step by step:
# 1. From F + L = 26 and L = 1.6F:
#    F + 1.6F = 26
#    2.6F = 26
#    F = 10
#    Therefore L = 16

F = 10
L = 16

# 2. From L = 3.2A:
#    16 = 3.2A
#    A = 5

A = 5

# 3. From A - E = 4:
E = A - 4
E = 1

# 4. From H - I = -4 and H + M = 26:
#    Let's try H = 16 (one of our available numbers)
H = 16
I = H + 4
I = 20

# 5. From I = 4J:
J = I / 4
J = 5

# 6. From D = 1.5H:
D = int(1.5 * H)
D = 24

# 7. From D - B = -14:
B = D + 14
B = 36

# 8. From H + M = 26:
M = 26 - H
M = 10

# 9. Remaining numbers: 2, 3, 7
# For C, G, K we can use these numbers
C = 7
G = 3
K = 2

solution = [A, B, C, D, E, F, G, H, I, J, K, L, M]
print(solution)