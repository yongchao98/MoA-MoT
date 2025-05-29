from itertools import product

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
    count = 0
    for truth_values in product([True, False], repeat=7):
        T = sum(truth_values)
        F = 7 - T
        
        # Check each statement
        s1 = T >= 4
        s2 = F <= 4
        s3 = T == 2
        s4 = F == 7
        s5 = (s3 or s4) and not (s3 and s4)
        s6 = is_prime(T)
        s7 = is_composite(F)
        
        # Check if the truth values satisfy the conditions
        if (s1 == truth_values[0] and
            s2 == truth_values[1] and
            s3 == truth_values[2] and
            s4 == truth_values[3] and
            s5 == truth_values[4] and
            s6 == truth_values[5] and
            s7 == truth_values[6]):
            count += 1
    return count

print(count_valid_solutions())