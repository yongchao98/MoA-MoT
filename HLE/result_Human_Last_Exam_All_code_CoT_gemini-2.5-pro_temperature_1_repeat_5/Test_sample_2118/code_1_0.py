# A dictionary to store computed coefficients for memoization
memo = {}

def get_coefficient(n):
    """
    Calculates the coefficient a_n of the x^-n term in the asymptotic expansion
    of f(x) using the derived recurrence relation.
    Uses memoization to avoid redundant calculations.
    """
    if n in memo:
        return memo[n]
    
    # Base case: The series starts at n=2. For the purpose of the recurrence,
    # we can define coefficients for n < 2 as 0.
    if n < 2:
        return 0
    
    # Recurrence relations derived from the functional equation:
    if n % 2 != 0:
        # For odd n (n >= 3), a_n = 1
        result = 1
    else: # n is even
        # For even n (n >= 2), a_n = 1 - a_{n/2}
        result = 1 - get_coefficient(n // 2)
        
    memo[n] = result
    return result

# We need to find the number of non-zero terms for n from 2 up to 100.
limit = 100
nonzero_count = 0
zero_count = 0

# Iterate from n=2 to 100 to check each term.
for n in range(2, limit + 1):
    if get_coefficient(n) != 0:
        nonzero_count += 1
    else:
        zero_count += 1

total_terms = limit - 1 

print("Calculation of the number of nonzero terms:")
print(f"Total terms considered (n=2 to 100): {total_terms}")
print(f"Number of terms with zero coefficients: {zero_count}")
print(f"Number of terms with nonzero coefficients: {nonzero_count}")
print("\nFinal Equation:")
print(f"{total_terms} - {zero_count} = {nonzero_count}")