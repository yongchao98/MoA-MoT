from itertools import permutations

def check_conditions(values):
    A, B, C, D, E, F, G, H, I, J, K = values
    
    # Check all conditions
    conditions = [
        H > G,                  # 1
        H + K == 273,          # 2
        B - F == 60,           # 3
        K - B == 145,          # 4
        F - I == 11,           # 5
        D - E == -59,          # 6
        B == 4.0 * F,          # 7
        G + H == 51,           # 8
        K == 3.0 * E,          # 9
        C + K == 264,          # 10
        I == 3.0 * G           # 11
    ]
    
    return all(conditions)

numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]
found = False

for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        found = True
        break

if not found:
    print("No solution found")