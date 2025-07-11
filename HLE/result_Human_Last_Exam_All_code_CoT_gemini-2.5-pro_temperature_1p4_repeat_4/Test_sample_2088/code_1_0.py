# The problem is to compute (12)^4 * (Integral)^4.
# We found the value of the integral to be 5^(1/4) / 12.
# Now we substitute this value back into the original expression.
# The calculation becomes (12)^4 * ( (5^(1/4) / 12) )^4.
# This simplifies to (12^4) * ( 5 / 12^4 ), which is 5.

val_12_to_4 = 12**4
numerator = 5
denominator = 12**4

# We are computing the expression: (12^4) * (5 / 12^4)
result = val_12_to_4 * (numerator / denominator)

print(f"The simplified expression to compute is ({12}^{4}) * ({5} / {12}^{4})")
print(f"= {val_12_to_4} * ({numerator} / {denominator})")
print(f"= {int(result)}")