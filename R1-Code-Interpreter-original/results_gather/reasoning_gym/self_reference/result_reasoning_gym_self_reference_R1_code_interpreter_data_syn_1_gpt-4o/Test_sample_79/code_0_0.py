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

def check_statements():
    possible_solutions = 0
    for true_count in range(8):
        false_count = 7 - true_count
        s1 = true_count >= 6
        s2 = false_count <= 1
        s3 = true_count == 2
        s4 = false_count == 3
        s5 = (s3 or s4) and not (s3 and s4)
        s6 = is_prime(true_count)
        s7 = is_composite(false_count)
        
        if (s1 + s2 + s3 + s4 + s5 + s6 + s7) == true_count:
            possible_solutions += 1
    return possible_solutions

print(check_statements())