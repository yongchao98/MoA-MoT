# Let's try to find a combination of operations to reach 147
numbers = [94, 90, 80, 81, 81, 72]

# Try different combinations
for a in numbers:
    for b in numbers:
        if a != b:
            for c in numbers:
                if c != a and c != b:
                    for d in numbers:
                        if d != a and d != b and d != c:
                            for e in numbers:
                                if e != a and e != b and e != c and e != d:
                                    for f in numbers:
                                        if f != a and f != b and f != c and f != d and f != e:
                                            # Check if any combination of operations results in 147
                                            if a + b - c + d - e + f == 147:
                                                print(f"{a} + {b} - {c} + {d} - {e} + {f}")
                                            elif a * b - c + d - e + f == 147:
                                                print(f"{a} * {b} - {c} + {d} - {e} + {f}")
                                            elif a + b * c - d + e - f == 147:
                                                print(f"{a} + {b} * {c} - {d} + {e} - {f}")
                                            elif a + b - c * d + e - f == 147:
                                                print(f"{a} + {b} - {c} * {d} + {e} - {f}")
                                            elif a + b - c + d * e - f == 147:
                                                print(f"{a} + {b} - {c} + {d} * {e} - {f}")
                                            elif a + b - c + d - e * f == 147:
                                                print(f"{a} + {b} - {c} + {d} - {e} * {f}")