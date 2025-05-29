def verify_solution():
    # The solution we got
    values = [96, 24, 16, 45, 7, 2, 3, 10, 28, 5, 36, 15, 50]
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Verify all constraints and print them
    print(f"1. D = 3.0L: {D} = 3.0 * {L} -> {D == 3.0 * L}")
    print(f"2. F + L = 17: {F} + {L} = {F + L} -> {F + L == 17}")
    print(f"3. I + J = 33: {I} + {J} = {I + J} -> {I + J == 33}")
    print(f"4. I - D = -17: {I} - {D} = {I - D} -> {I - D == -17}")
    print(f"5. I - E = 21: {I} - {E} = {I - E} -> {I - E == 21}")
    print(f"6. D - C = 29: {D} - {C} = {D - C} -> {D - C == 29}")
    print(f"7. A = 4.0B: {A} = 4.0 * {B} -> {A == 4.0 * B}")
    print(f"8. F - C = -14: {F} - {C} = {F - C} -> {F - C == -14}")
    print(f"9. L = 3.0J: {L} = 3.0 * {J} -> {L == 3.0 * J}")
    print(f"10. L - D = -30: {L} - {D} = {L - D} -> {L - D == -30}")
    print(f"11. F - B = -22: {F} - {B} = {F - B} -> {F - B == -22}")
    
    # Print sorted list of used numbers to verify all numbers are used
    print("\nUsed numbers (sorted):", sorted(values))
    print("Should be: [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]")

verify_solution()