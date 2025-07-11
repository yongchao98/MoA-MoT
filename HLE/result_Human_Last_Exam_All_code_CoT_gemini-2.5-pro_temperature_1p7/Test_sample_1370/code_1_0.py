# The problem asks for the maximum number of mutually independent events (m)
# for an experiment of rolling 100 regular 6-sided dice.

# 1. Determine the size of the sample space.
# For 100 dice with 6 sides each, the total number of outcomes is 6^100.
num_dice = 100
sides = 6
# Sample space size |Ω| = sides ^ num_dice = 6^100.

# 2. Use the theorem on the maximum number of independent events.
# For a finite equiprobable sample space of size N, the maximum number of
# mutually independent events is the sum of the exponents in the prime factorization of N.

# 3. Find the prime factorization of the sample space size.
# N = 6^100 = (2 * 3)^100 = 2^100 * 3^100.

# 4. Identify the exponents of the prime factors.
# The prime factors are 2 and 3.
# The exponent for the prime factor 2 is 100.
exponent_for_2 = 100
# The exponent for the prime factor 3 is 100.
exponent_for_3 = 100

# 5. Sum the exponents to get the maximum number of events, m.
m = exponent_for_2 + exponent_for_3

# 6. Print the final calculation and the result.
# The problem asks to output each number in the final equation.
print(f"{exponent_for_2} + {exponent_for_3} = {m}")