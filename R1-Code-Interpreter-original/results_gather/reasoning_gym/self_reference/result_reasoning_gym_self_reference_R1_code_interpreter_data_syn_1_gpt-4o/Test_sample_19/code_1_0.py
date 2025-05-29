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

solutions = 0

for true_count in range(8):
    false_count = 7 - true_count
    
    # Check each statement
    s1 = (true_count >= 7)
    s2 = (true_count >= 6)
    s3 = (true_count == 1)
    s4 = (false_count == 7)
    s5 = (s3 != s4)
    s6 = is_prime(true_count)
    s7 = is_composite(false_count)
    
    # Count valid solutions
    if [s1, s2, s3, s4, s5, s6, s7].count(True) == true_count:
        solutions += 1

print(solutions)