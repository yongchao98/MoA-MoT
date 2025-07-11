import math

def solve_and_print():
    """
    This function determines the number of nonzero terms up to and including x^-100
    in the asymptotic expansion of the function f(x) defined by the given equation.

    The condition for a coefficient c_n to be non-zero is determined by writing n
    as n = 2^p * m, where m is an odd number.
    - If m = 1 (n is a power of 2), c_n is non-zero if and only if p is odd.
    - If m > 1, c_n is non-zero if and only if p is even.

    The code counts the number of integers n from 2 to 100 that satisfy these
    conditions and prints the breakdown of the count.
    """
    limit = 100
    
    # List to store the counts for each case
    counts = []

    # Case 1: n = 2^p where p is odd
    count_m1_p_odd = 0
    p = 1
    while True:
        n = 2**p
        if n > limit:
            break
        # Check if p is odd
        if p % 2 != 0:
            count_m1_p_odd += 1
        p += 1
    counts.append(count_m1_p_odd)

    # Case 2: n = 2^p * m where m > 1 is odd and p is even
    max_p = int(math.log2(limit))
    for p in range(max_p + 1):
        # Check if p is even
        if p % 2 == 0:
            count_p_even = 0
            base = 2**p
            # Start with the smallest odd m > 1, which is 3
            m = 3
            while True:
                n = base * m
                if n > limit:
                    break
                count_p_even += 1
                # Move to the next odd number
                m += 2
            
            if count_p_even > 0:
                counts.append(count_p_even)

    total_count = sum(counts)
    
    # Format the counts into a final equation string
    equation_str = " + ".join(map(str, counts))
    
    print(f"The total number of nonzero terms is {total_count}.")
    print("This is calculated by summing the counts from the different cases:")
    print(f"{total_count} = {equation_str}")

solve_and_print()