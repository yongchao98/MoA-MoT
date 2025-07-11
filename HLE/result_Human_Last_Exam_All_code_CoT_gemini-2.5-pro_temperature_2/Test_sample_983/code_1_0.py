def get_sum_of_distinct_prime_factors(n):
    """
    Calculates the sum of distinct prime factors of a given integer.
    """
    factors = set()
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.add(temp)
    return sum(factors)

# The last number in the sequence
last_number = 2352

# Calculate the sum of its distinct prime factors
difference = get_sum_of_distinct_prime_factors(last_number)

# Calculate the next number in the sequence
next_number = last_number + difference

# Print the final equation
# To show the list of distinct prime factors for clarity in the output
distinct_factors = set()
d = 2
temp = last_number
while d * d <= temp:
    if temp % d == 0:
        distinct_factors.add(d)
        while temp % d == 0:
            temp //= d
    d += 1
if temp > 1:
    distinct_factors.add(temp)

print(f"The next number is calculated by adding the sum of the distinct prime factors of the last number to itself.")
print(f"The distinct prime factors of {last_number} are {sorted(list(distinct_factors))}.")
print(f"Their sum is {' + '.join(map(str, sorted(list(distinct_factors))))} = {difference}.")
print(f"So, the next number is {last_number} + {difference} = {next_number}.")