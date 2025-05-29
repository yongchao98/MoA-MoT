def verify_solution():
    # The solution we got
    A, B, C, D, E, F, G, H, I, J, K = 75, 9, 121, 3.0, 48.0, 36.0, 27, 150, 225.0, 16, 80.0
    
    # Check all constraints
    print("1. I = 1.5H:", 225.0 == 1.5 * 150)
    print("2. F = 4.0B:", 36.0 == 4.0 * 9)
    print("3. C + E = 169:", 121 + 48.0 == 169)
    print("4. E = 3.0J:", 48.0 == 3.0 * 16)
    print("5. J - G = -11:", 16 - 27 == -11)
    print("6. C > F:", 121 > 36.0)
    print("7. I + K = 305:", 225.0 + 80.0 == 305)
    print("8. I - D = 222:", 225.0 - 3.0 == 222)
    print("9. J - F = -20:", 16 - 36.0 == -20)
    print("10. F + H = 186:", 36.0 + 150 == 186)
    
    # Check if all numbers are from the given set
    given_set = set([3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225])
    solution_set = set([75, 9, 121, 3.0, 48.0, 36.0, 27, 150, 225.0, 16, 80.0])
    print("\nAll numbers in given set:", all(int(x) in given_set for x in solution_set))

print("Solution verification:")
verify_solution()