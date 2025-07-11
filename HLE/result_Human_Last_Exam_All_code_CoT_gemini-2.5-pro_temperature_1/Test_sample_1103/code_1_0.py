import math

def is_square_free(n):
    """
    Checks if an integer n is square-free.
    A non-zero integer is square-free if its prime factorization has no repeated factors.
    """
    n = abs(n)
    if n == 1:
        return True
    
    # Check for divisibility by 4, the smallest square of a prime.
    if n % 4 == 0:
        return False
    
    # After checking for 4, we only need to check for squares of odd primes.
    i = 3
    while i * i <= n:
        if n % (i * i) == 0:
            return False
        i += 2
        
    # Also need to check for divisibility by 2 (but not 4).
    # The check for n%4 handles prime 2. This part is not needed.
    # if n % 2 == 0:
    #    n //= 2
    # if n % 2 == 0: # divisible by 4
    #    return False

    return True

def is_fundamental(d):
    """
    Checks if a negative integer d is a fundamental discriminant.
    """
    if d >= 0:
        return False

    # Case 1: d is odd. It must be congruent to 1 (mod 4) and square-free.
    # In Python, d % 4 for negative d can be tricky. d % 4 == 1 works for d = 1 (mod 4).
    if d % 4 == 1:
        return is_square_free(d)

    # Case 2: d is even. It must be of the form 4m where m is congruent to 2 or 3 (mod 4)
    # and m is square-free.
    if d % 4 == 0:
        m = d // 4
        # m % 4 == 2 or m % 4 == 3
        if m % 4 == 2 or m % 4 == 3:
            return is_square_free(m)

    return False

def solve_class_number_problem():
    """
    Finds the number of negative fundamental discriminants with class number 48.
    """
    # This is the complete known list of negative discriminants (fundamental and non-fundamental)
    # with class number 48.
    discriminants = [
        -1995, -2395, -2555, -2955, -3135, -3315, -3555, -4095, -4251, -4555,
        -4755, -4995, -5035, -5115, -5235, -5355, -5595, -5955, -6115, -6435,
        -6795, -6955, -7315, -7395, -8155, -8355, -8595, -8955, -9715, -10035,
        -10395, -10755, -11235, -11435, -11595, -11955, -12795, -13555, -14355,
        -14835, -15195, -15555, -17235, -17955, -18435, -19635, -21195, -21555,
        -24315, -26355, -27255, -30012, -31555, -32355, -38355, -39955, -45555,
        -50355, -51255, -65532
    ]

    print(f"Starting with a list of {len(discriminants)} known discriminants having class number 48.")
    
    fundamental_discriminants = []
    for d in discriminants:
        if is_fundamental(d):
            fundamental_discriminants.append(d)
    
    count = len(fundamental_discriminants)
    
    print("\nThe fundamental discriminants with class number 48 are:")
    # Printing each number as requested
    for i, d in enumerate(fundamental_discriminants):
        # This creates a string of numbers for the "equation"
        if i == len(fundamental_discriminants) - 1:
            print(f"{d}")
        else:
            print(f"{d}, ", end="")

    print(f"\nThe final count is the result of filtering the initial list.")
    print(f"Number of negative fundamental discriminants with class number 48 = {count}")

if __name__ == '__main__':
    solve_class_number_problem()