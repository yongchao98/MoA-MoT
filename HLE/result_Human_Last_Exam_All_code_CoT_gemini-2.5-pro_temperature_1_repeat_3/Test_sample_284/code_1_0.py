# The moment curve in R^3 is given by gamma(t) = (t^d1, t^d2, t^d3)
# In this problem, the exponents are d1=1, d2=2, d3=3.
d1 = 1
d2 = 2
d3 = 3

# The critical exponent p is the sum of these exponents.
# This comes from the scaling argument where the volume of the
# region of constructive interference scales as delta to the power of
# -(d1 + d2 + d3).
p_critical = d1 + d2 + d3

print(f"The critical exponent is given by the sum of the powers of the moment curve.")
print(f"The equation is: {d1} + {d2} + {d3} = {p_critical}")
print(f"Thus, the largest possible value of p is {p_critical}.")
