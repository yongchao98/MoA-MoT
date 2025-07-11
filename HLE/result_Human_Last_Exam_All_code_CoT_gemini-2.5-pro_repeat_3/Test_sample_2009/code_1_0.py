import math

# We found the squared length of AE to be 578/15
# AE^2 = 578 / 15
# We need to express AE in the form a*sqrt(b)/c

# We can find the integer part of the square root of 578
# 24^2 = 576
# So sqrt(578) = sqrt(2 * 289) = 17 * sqrt(2)
# AE = (17 * sqrt(2)) / sqrt(15)
# To rationalize the denominator, multiply numerator and denominator by sqrt(15)
# AE = (17 * sqrt(2) * sqrt(15)) / 15
# AE = (17 * sqrt(30)) / 15

# So a=17, b=30, c=15
a = 17
b = 30
c = 15

# Check conditions: b is square-free.
# Prime factors of 30 are 2, 3, 5. No squared factors. So it is square-free.
# Check gcd(a,c) = 1
# gcd(17, 15) = 1, since 17 is a prime and not a factor of 15.

# The problem asks for a+b+c
result = a + b + c

print(f"The length of AE is represented as (a*sqrt(b))/c.")
print(f"From the derivation, we have a = {a}, b = {b}, c = {c}.")
print(f"The equation for AE is: AE = ({a}*sqrt({b}))/{c}")
print(f"The sum a+b+c is {a} + {b} + {c} = {result}")