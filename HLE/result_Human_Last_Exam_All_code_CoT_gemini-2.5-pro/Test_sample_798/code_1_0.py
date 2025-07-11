import math

def is_prime(n):
    """
    Checks if a number is prime.
    """
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

def solve_partition_problem():
    """
    Calculates the minimum number of distinct-distance sets needed to partition
    the integers from 10001 to 42149572.
    """
    start = 10001
    end = 42149572

    # Step 1 & 2: Calculate the total number of integers, N.
    N = end - start + 1
    
    print(f"The range of integers is from {start} to {end}.")
    print(f"The total number of integers to partition is N = {end} - {start} + 1 = {N}.")
    print("-" * 30)

    # Step 3: Find q such that N = q * (q + 1).
    # This is equivalent to solving the quadratic equation q^2 + q - N = 0.
    # The positive root is q = (-1 + sqrt(1 + 4*N)) / 2.
    delta = 1 + 4 * N
    q_float = (-1 + math.sqrt(delta)) / 2
    q = int(q_float)

    # Verify that q is an integer and that N = q * (q + 1).
    if q_float == q and N == q * (q + 1):
        print(f"We observe that N has a special form: N = q * (q + 1).")
        print(f"Solving q^2 + q - {N} = 0 for q, we find q = {q}.")
        print(f"So, the equation is {N} = {q} * ({q} + 1) = {q} * {q+1}.")
        print("-" * 30)
    else:
        print("N does not have the expected form q*(q+1). The method is not applicable.")
        return

    # Step 4 & 5: Apply number theory results.
    # Check if q is a prime number (a prime is a prime power).
    if is_prime(q):
        print(f"The value q = {q} is a prime number.")
        print("For N = q*(q+1) where q is a prime power, it is known that:")
        print(f"1. The largest possible distinct-distance set in {{1, ..., {N}}} has size q+1 = {q+1}.")
        print(f"2. The minimum number of sets 'k' must be at least N / (q+1) = {N} / {q+1} = {q}.")
        print(f"3. It is also known that {{1, ..., {N}}} can be partitioned into exactly q = {q} distinct-distance sets.")
        print("\nCombining these facts, the minimum number of sets required is q.")
        
        final_answer = q
        print(f"\nThe final answer is {final_answer}.")
    else:
        print(f"The value q = {q} is not a prime power. This specific method does not apply.")

solve_partition_problem()
<<<6491>>>