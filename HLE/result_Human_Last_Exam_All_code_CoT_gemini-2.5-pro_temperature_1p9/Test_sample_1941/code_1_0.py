import collections
import math

def is_prime(n):
    """
    Checks if a number n is prime using an optimized trial division method.
    A number is prime if it is greater than 1 and has no positive divisors other than 1 and itself.
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

# --- Main execution ---
# Step 1: Find the first 1,000 prime numbers with only digits 0 and 1 (Set C).
set_c = []
# Use a deque as an efficient queue for breadth-first generation of '01' numbers.
# We start with '1' as the seed for all other numbers.
q = collections.deque(['1'])

target_count = 1000

print("Generating the first 1,000 primes containing only digits 0 and 1...")

while len(set_c) < target_count:
    # Get the next candidate number string from the queue.
    s_num = q.popleft()

    # Convert to integer and test for primality.
    # The first number we check from the queue is '1', which is not prime.
    # The first actual prime found will be 11.
    num = int(s_num)
    if is_prime(num):
        set_c.append(s_num)

    # Enqueue the next potential candidates by appending '0' and '1'.
    # This ensures we are generating and testing numbers in increasing order.
    q.append(s_num + '0')
    q.append(s_num + '1')

# Step 2: Analyze the found primes based on the Game of Life rule.
# A one-cell-thick pattern dies out if and only if its decimal representation
# does not contain the substring "111".
die_out_count = 0
for p_str in set_c:
    if "111" not in p_str:
        die_out_count += 1

# Step 3: Output the result.
# The 'final equation' is the summary of our findings.
total_primes_in_c = len(set_c)

print(f"\nThe total number of primes found in set C is: {total_primes_in_c}")
print(f"The number of these primes that represent patterns which will ultimately die out is: {die_out_count}")
print("<<<" + str(die_out_count) + ">>>")