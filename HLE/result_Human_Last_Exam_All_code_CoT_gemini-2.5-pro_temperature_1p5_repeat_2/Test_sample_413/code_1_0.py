from fractions import Fraction

# Step 1 & 2: Define coefficients a and b based on m and alpha
# We derived that m must be an integer >= 0. Let's take m=0.
m = 0
# We derived alpha = -(4m+5)/8
alpha = -(4 * m + 5) / Fraction(8)

# The roots are m, m+1, and alpha
root1 = m
root2 = m + 1
root3 = alpha

# Step 3: Define f(x) and calculate the coefficients a,b,c
# f(x) = (x-r1)(x-r2)(x-r3)
# f(x) = x^3 - (r1+r2+r3)x^2 + (r1r2+r1r3+r2r3)x - r1r2r3
a = -(root1 + root2 + root3)
b = root1*root2 + root1*root3 + root2*root3
c = -(root1 * root2 * root3)

# We can now write f(x) explicitly
# f(x) = x^3 + a*x^2 + b*x + c

# Step 4: Compute f(3)
x = 3
f_3 = x**3 + a * x**2 + b * x + c

# The final calculation is f(3) = 27 + a*9 + b*3 + c
term1 = 27
term2 = a * 9
term3 = b * 3
term4 = c

# Output the components of the final sum and the result
# The problem asks for the exact value of f(3).
# Our f(x) is x^3 - 3/8*x^2 - 5/8*x
# So f(3) = 3^3 - (3/8)*3^2 - (5/8)*3
val1 = 3**3
val2 = Fraction(3,8) * 3**2
val3 = Fraction(5,8) * 3

print(f"The function is f(x) = x^3 + ({a})x^2 + ({b})x + ({c})")
print(f"To compute f(3), we evaluate 3^3 + ({a})*3^2 + ({b})*3 + ({c})")
print(f"f(3) = {val1} - {val2} - {val3}")
print(f"f(3) = {val1} - ({val2+val3})")
print(f"f(3) = {val1 - (val2+val3)}")
print(f"The exact value of f(3) is {f_3}")
