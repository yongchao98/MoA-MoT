import collections

def solve():
    """
    This function solves the problem by finding the first 1,000 primes
    made of only 0s and 1s, and then counting how many of them would
    die out in Conway's Game of Life.
    """

    def is_prime(n):
        """
        Checks if a number n is prime using an optimized trial division method.
        """
        if n < 2:
            return False
        if n == 2 or n == 3:
            return True
        if n % 2 == 0 or n % 3 == 0:
            return False
        i = 5
        while i * i <= n:
            if n % i == 0 or n % (i + 2) == 0:
                return False
            i += 6
        return True

    # Step 1 & 2: Generate candidate numbers and test for primality to find set C.
    C = []
    # Use a deque for an efficient queue in our breadth-first search.
    queue = collections.deque(["1"])

    # The first prime with only 0s and 1s is 11.
    # The generator will check 1 (not prime), then 10, then 11.
    while len(C) < 1000:
        s = queue.popleft()
        n = int(s)

        # Check for primality. n=1 is not prime.
        if n > 1 and is_prime(n):
            C.append(n)

        # Continue the generation by appending '0' and '1'.
        queue.append(s + "0")
        queue.append(s + "1")

    # Step 3: Count how many members of C will die out.
    # A pattern dies out if its string representation does not contain "111".
    die_out_count = 0
    for p in C:
        if "111" not in str(p):
            die_out_count += 1
    
    # The "final equation" is the count of numbers that meet the criteria.
    # We print each part of the "equation", which in this case is just the final number.
    print(die_out_count)

solve()