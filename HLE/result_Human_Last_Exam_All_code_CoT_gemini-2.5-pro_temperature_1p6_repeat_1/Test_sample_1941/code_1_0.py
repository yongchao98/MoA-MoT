import math

def is_prime(n):
    """Checks if a number is prime."""
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True

def has_only_01_digits(n):
    """Checks if a number's decimal representation contains only '0' and '1'."""
    return all(c in '01' for c in str(n))

def find_dying_patterns():
    """
    Finds the first 1,000 primes with only digits 0 and 1,
    and counts how many will die out in Conway's Game of Life.
    """
    set_c = []
    num = 2
    # Step 1: Generate Set C, the first 1,000 primes with only '0' and '1' digits.
    while len(set_c) < 1000:
        if has_only_01_digits(num) and is_prime(num):
            set_c.append(num)
        num += 1

    died_out_count = 0
    patterns_that_die = []

    # Step 2: Analyze each number in Set C based on the "die out" condition.
    # A pattern dies out if its decimal representation does not contain "111".
    for prime in set_c:
        if '111' not in str(prime):
            died_out_count += 1
            patterns_that_die.append(prime)

    # Step 3: Print the final result.
    print("This program calculates how many of the first 1,000 primes containing only digits 0 and 1")
    print("would die out if encoded as one-cell-thick patterns in Conway's Game of Life.")
    print("\nA pattern dies out if its number representation does not contain the substring '111'.")
    print(f"\nOut of the {len(set_c)} primes in the set, a total of {died_out_count} patterns will die out.")
    print(f"\nThe equation is the count of primes found without the '111' substring:")
    print(f"{len(patterns_that_die)} = {died_out_count}")


if __name__ == '__main__':
    find_dying_patterns()
<<<831>>>