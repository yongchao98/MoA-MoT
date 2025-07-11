# The problem asks for the exponent alpha in the upper bound N^alpha for the measure of the set X.
# Based on the analytical derivation, the exponent alpha is given by the expression:
# alpha = 5/2 - 3 * (3/8)

a_num = 5
a_den = 2
b_num = 3
b_den = 8

# The expression for alpha simplifies to 5/2 - 9/8.
# To perform the subtraction, we find a common denominator, which is 8.
c_num = 9
c_den = 8

final_numerator = a_num * (c_den // a_den) - c_num
final_denominator = c_den

print(f"The formula for the exponent alpha is alpha = {a_num}/{a_den} - 3 * ({b_num}/{b_den}).")
print(f"This simplifies to the equation: alpha = {a_num}/{a_den} - {c_num}/{c_den}.")
print(f"Expressing with a common denominator gives: alpha = ({a_num * (c_den // a_den)})/{c_den} - {c_num}/{c_den}.")
print(f"The final calculation is: alpha = ({a_num * (c_den // a_den) - c_num})/{c_den} = {final_numerator}/{final_denominator}.")

alpha_value = final_numerator / final_denominator
print(f"The numerical value of alpha is {alpha_value}.")