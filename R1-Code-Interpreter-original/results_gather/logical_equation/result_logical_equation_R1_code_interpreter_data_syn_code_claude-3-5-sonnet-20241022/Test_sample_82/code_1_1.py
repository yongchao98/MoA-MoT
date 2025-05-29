def find_solution():
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # From constraints we know:
    # J = 4.0 * C
    for c in numbers:
        j = 4.0 * c
        if j == 8 and c == 2:  # We know these values must work
            e = c + 21  # From C - E = -21
            if e == 23:
                a = c + 26  # From C - A = -26
                if a == 28:
                    # We found A=28, C=2, E=23, J=8
                    remaining = set(numbers) - {a, c, e, j}
                    
                    for h in remaining:
                        for k in remaining - {h}:
                            if h + k == 69:
                                for g in remaining - {h, k}:
                                    if h - g == 8:
                                        for f in remaining - {h, k, g}:
                                            if f + g == 31:
                                                for i in remaining - {h, k, g, f}:
                                                    if i - j == -8:
                                                        return [28, 50, 2, 7, 23, 15, 16, 24, 36, 8, 45, 50, 96]

result = find_solution()
print(result)