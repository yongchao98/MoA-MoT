def get_prime_factorization(num):
    """
    Returns the prime factorization of a number as a dictionary of primes to their exponents.
    """
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def is_sum_of_two_squares(n):
    """
    Checks if an integer n can be written as the sum of two squares.
    This is true if all prime factors of the form 4k+3 have an even exponent.
    """
    if n == 0 or n == 1:
        return True
    
    factors = get_prime_factorization(n)
    for p, exp in factors.items():
        if p % 4 == 3:
            if exp % 2 != 0:
                return False
    return True

def find_subset():
    """
    Finds the subset of integers for which the number of n x n tilings is always even.
    The condition is: t=2 or t is not a sum of two squares.
    """
    input_set = [2, 3, 4, 5, 7, 9, 15]
    solution_subset = []
    
    print("Analyzing the set T = {2, 3, 4, 5, 7, 9, 15}\n")
    
    for t in input_set:
        is_sots = is_sum_of_two_squares(t)
        
        # The property holds for t=2 and for t which are NOT sums of two squares.
        if t == 2 or not is_sots:
            solution_subset.append(t)
            print(f"t = {t}: Is a sum of two squares? {is_sots}. It BELONGS to the subset.")
        else:
            print(f"t = {t}: Is a sum of two squares? {is_sots}. It DOES NOT belong to the subset.")

    print("\nThe subset of integers for which the statement is true is:")
    print(solution_subset)

if __name__ == '__main__':
    find_subset()