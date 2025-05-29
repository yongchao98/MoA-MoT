def find_solution():
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # Find possible M,K pairs where K = 4.5M
    for M in numbers:
        K = 4.5 * M
        if K in numbers:
            # Find possible G,E pairs where E = 1.5G
            for G in numbers:
                E = 1.5 * G
                if E in numbers and E + K == 81:  # Check E + K = 81
                    # Find possible C,F pairs where F = 3.2C
                    for C in numbers:
                        F = 3.2 * C
                        if F in numbers and C - E == -31:  # Check C - E = -31
                            # Find possible B,J pairs where B + J = 5
                            for B in numbers:
                                for J in numbers:
                                    if B + J == 5 and K > J and M > B:
                                        L = 1.4 * C
                                        if L in numbers:
                                            remaining = set(numbers) - {B,C,E,F,G,J,K,L,M}
                                            for I in remaining:
                                                if F + I == 112:
                                                    remaining2 = remaining - {I}
                                                    for A in remaining2:
                                                        if A + M == 38:
                                                            remaining3 = remaining2 - {A}
                                                            D = remaining3.pop()
                                                            result = [A,B,C,D,E,F,G,H,I,J,K,L,M]
                                                            if len(set(result)) == 13:  # Check all numbers are different
                                                                print([A,B,C,D,E,F,G,H,I,J,K,L,M])
                                                                return

find_solution()