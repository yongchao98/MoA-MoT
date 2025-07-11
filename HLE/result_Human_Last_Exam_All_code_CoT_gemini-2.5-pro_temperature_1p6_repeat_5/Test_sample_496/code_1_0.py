# This script calculates the total rank of the equivariant cohomology ring A
# in degrees up to 100.
# The Poincare series for the rank of A is P(t) = (t^4 + t^7) / (1-t^4)^2.
# We expand this as a power series and sum the coefficients for degrees <= 100.
# P(t) = sum_{k=0 to inf} (k+1)*t^(4k+4) + sum_{k=0 to inf} (k+1)*t^(4k+7)

total_rank = 0

# Calculate the contribution from the first series: sum_{k=0 to inf} (k+1)*t^(4k+4)
# We need degrees 4k+4 <= 100, which means 4k <= 96, so k <= 24.
sum1 = 0
for k in range(25):  # k from 0 to 24
    rank = k + 1
    sum1 += rank

print(f"Sum of ranks from the first series (degrees 4k+4): {sum1}")

# Calculate the contribution from the second series: sum_{k=0 to inf} (k+1)*t^(4k+7)
# We need degrees 4k+7 <= 100, which means 4k <= 93, so k <= 23.
sum2 = 0
for k in range(24):  # k from 0 to 23
    rank = k + 1
    sum2 += rank

print(f"Sum of ranks from the second series (degrees 4k+7): {sum2}")

total_rank = sum1 + sum2
print(f"Total rank up to degree 100 is: {sum1} + {sum2} = {total_rank}")

# Final Answer
# print(total_rank)