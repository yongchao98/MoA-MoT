def is_prime(n):
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

possible_solutions = 0

for true_count in range(8):
    false_count = 7 - true_count
    
    # Check each condition
    if true_count >= 1 and false_count <= 3:
        if (true_count == 4) != (false_count == 3):
            if is_prime(true_count) and is_composite(false_count):
                possible_solutions += 1

print(possible_solutions)