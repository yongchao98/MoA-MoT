import sympy

# Set up the symbolic variable
t = sympy.Symbol('t', real=True)

# Define the y-component of the curve, y = |t^3|
y = sympy.Abs(t**3)

print("Analyzing the smoothness of y(t) = |t^3| at t=0.")

# Calculate the first derivative
# For t>0, y=t^3, y'=3t^2. For t<0, y=-t^3, y'=-3t^2. At t=0, y'(0)=0.
# The expression is 3*t*Abs(t).
y_prime = sympy.diff(y, t)
print(f"The first derivative y'(t) is: {y_prime}")
y_prime_val_at_0 = sympy.limit(y_prime, t, 0)
print(f"The value of y'(t) at t=0 is: {y_prime_val_at_0}")


# Calculate the second derivative
# For t>0, y''=6t. For t<0, y''=-6t. At t=0, y''(0)=0.
# The expression is 6*Abs(t).
y_prime_prime = sympy.diff(y_prime, t)
print(f"The second derivative y''(t) is: {y_prime_prime}")
y_prime_prime_val_at_0 = sympy.limit(y_prime_prime, t, 0)
print(f"The value of y''(t) at t=0 is: {y_prime_prime_val_at_0}")


# Calculate the third derivative
# For t>0, y'''=6. For t<0, y'''=-6.
# The derivative does not exist at t=0.
y_ppp_right = sympy.limit(sympy.diff(6*t, t), t, 0, dir='+')
y_ppp_left = sympy.limit(sympy.diff(-6*t, t), t, 0, dir='-')

print(f"The limit of the third derivative as t->0 from the right is: {y_ppp_right}")
print(f"The limit of the third derivative as t->0 from the left is: {y_ppp_left}")

print("\nSince the left and right limits of the third derivative at t=0 are not equal, the third derivative does not exist.")
print("This demonstrates that the curve gamma(t) = (t^3, |t^3|) is not a C-infinity (smooth) curve, which supports the conclusion that statement B is false.")
print("However, statement C is false for more fundamental reasons regarding dimension and topology.")
