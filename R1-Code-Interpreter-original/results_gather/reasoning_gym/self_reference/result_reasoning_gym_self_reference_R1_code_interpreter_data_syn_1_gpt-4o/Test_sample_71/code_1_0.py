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

def count_valid_solutions():
    valid_solutions = 0
    for true_count in range(8):
        false_count = 7 - true_count
        statement1 = true_count >= 4
        statement2 = false_count <= 4
        statement3 = true_count == 4
        statement4 = false_count == 3
        statement5 = (statement3 or statement4) and not (statement3 and statement4)
        statement6 = is_prime(true_count)
        statement7 = is_composite(false_count)
        
        if (statement1 and statement2 and statement5 and statement6 and statement7):
            valid_solutions += 1
    return valid_solutions

print(count_valid_solutions())