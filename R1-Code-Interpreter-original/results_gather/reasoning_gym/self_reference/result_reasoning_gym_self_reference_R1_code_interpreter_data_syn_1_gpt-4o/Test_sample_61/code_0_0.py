def is_prime(n):
    if n <= 1:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

possible_solutions = 0

for T in range(8):  # T can be from 0 to 7
    F = 7 - T
    if (T >= 3 and  # Statement 1
        T == 3 and  # Statement 2 is false
        T != 4 and  # Statement 3 is false
        T != 4 and  # Statement 4 is false
        (T == 3 or T == 4) and not (T == 3 and T == 4) and  # Statement 5
        is_prime(T) and  # Statement 6
        is_composite(F)):  # Statement 7
        possible_solutions += 1

print(possible_solutions)