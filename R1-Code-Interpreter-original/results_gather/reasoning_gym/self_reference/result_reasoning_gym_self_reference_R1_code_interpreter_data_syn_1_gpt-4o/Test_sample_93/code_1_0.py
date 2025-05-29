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
    
    # Check each statement
    s1 = (true_count >= 7)
    s2 = (false_count <= 6)
    s3 = (true_count == 5)
    s4 = (false_count == 6)
    s5 = (s3 != s4)
    s6 = is_prime(true_count)
    s7 = is_composite(false_count)
    
    # Check if all conditions are satisfied
    if s5 and s6 and s7:
        possible_solutions += 1

print(possible_solutions)