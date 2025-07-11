import math
from fractions import Fraction

# We are looking for sets of 4 distinct positive integers {a, b, c, d}
# such that 1/a + 1/b + 1/c + 1/d = 1.
# T(4) is the sum of all possible values of a+b+c+d.

solutions = []
sums_set = set()

# We assume a < b < c < d to find unique sets.
# This implies 1 = 1/a + 1/b + 1/c + 1/d < 4/a, so a < 4.
# a=1 is not possible, so a can be 2 or 3.
# We iterate through the possible values.

# Loop for a
for a in range(2, 4):
    # Determine search range for b
    # 1 - 1/a = 1/b + 1/c + 1/d < 3/b  =>  (a-1)/a < 3/b  =>  b < 3a/(a-1)
    limit_b = math.ceil(3 * a / (a - 1))
    for b in range(a + 1, limit_b):
        sum_ab = Fraction(1, a) + Fraction(1, b)
        if sum_ab >= 1:
            continue
        
        # Determine search range for c
        # 1 - sum_ab = 1/c + 1/d < 2/c  => c < 2 / (1 - sum_ab)
        rem_ab = Fraction(1, 1) - sum_ab
        limit_c = math.ceil(2 / rem_ab)
        for c in range(b + 1, limit_c):
            sum_abc = sum_ab + Fraction(1, c)
            if sum_abc >= 1:
                continue
            
            # Calculate d from the remaining fraction
            rem_abc = Fraction(1, 1) - sum_abc
            
            # A valid solution for d exists if the remaining fraction is of the form 1/d
            if rem_abc.numerator == 1:
                d = rem_abc.denominator
                # The final check is that d must be distinct and larger than c
                if d > c:
                    solution_set = (a, b, c, d)
                    current_sum = a + b + c + d
                    solutions.append(solution_set)
                    sums_set.add(current_sum)

# Sort solutions for a clean and ordered output, based on their sum
solutions.sort(key=lambda s: sum(s))

print("The sets of four distinct positive integers (a, b, c, d) and their sums are:")
for s in solutions:
    a, b, c, d = s
    current_sum = sum(s)
    print(f"For the equation 1/{a} + 1/{b} + 1/{c} + 1/{d} = 1, the sum is {a} + {b} + {c} + {d} = {current_sum}")

total_sum = sum(sums_set)
print("\nThe set of unique sums, S(4), is:", sorted(list(sums_set)))
print(f"The sum of all elements in S(4), which is T(4), is: {total_sum}")
