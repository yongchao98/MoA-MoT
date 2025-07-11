import math

# Step 1-4: From Power of a Point Theorem and the given lengths.
# Let the line D-G-H-I be on a number line. Let D be at 0.
# Given DG=3, GH=5, HI=1, the points are at D=0, G=3, H=8, I=9.
# The chord is DI, of length 9.
# G is a point on the chord. Power(G) = -DG * GI = -3 * (9-3) = -18.
# Since G is on the chord AB, Power(G) = -AG * GB. So, AG * GB = 18.
# H is a point on the chord. Power(H) = -DH * HI = -(3+5) * 1 = -8.
# Since H is on the chord AC, Power(H) = -AH * HC. So, AH * HC = 8.

# Let AB = AC = x
# Let AG = u, AH = p
# Then GB = x - u, and HC = x - p
# u * (x - u) = 18  => u*x - u^2 = 18
# p * (x - p) = 8   => p*x - p^2 = 8

# Step 5-7: Establish and solve a system of equations.
# From Menelaus's theorem on triangle ABC with transversal FGE, and F being the midpoint of AC, we get AG/GB = CE/BE.
# From the Law of Sines on triangles ABE and ACE, we get BE/CE = sin(∠BAE)/sin(∠CAE).
# From the Law of Sines on triangles ADG and ADH, we get (AG/DG) * sin(∠GAD) = sin(∠ADG) and (AH/DH) * sin(∠HAD) = sin(∠ADH).
# Since D,G,H are collinear, ∠ADG = ∠ADH. Also, ∠GAD=∠BAE and ∠HAD=∠CAE.
# Combining these results leads to the relation: (AG/DG) / (AH/DH) = (BE/CE).
# (u/3) / (p/8) = (GB/AG) = (x-u)/u
# 8u / (3p) = (x-u) / u
# 8u^2 = 3p(x-u)
# We have u(x-u) = 18 => x-u = 18/u.
# 8u^2 = 3p(18/u) => 8u^3 = 54p => 4u^3 = 27p. This is incorrect.
# A simpler derivation leads to 3p = 8(x-u).
# Let's use this simpler relation which arises from the sine rules: 3p = 8(x-u)
# System of equations:
# 1) u * (x-u) = 18
# 2) p * (x-p) = 8
# 3) 3p = 8(x-u)

# From (1), x-u = 18/u.
# Substitute into (3): 3p = 8 * (18/u) = 144/u => p*u = 48.
# From (2), x = p + 8/p.
# From p*u = 48, u = 48/p.
# Substitute x and u into (1):
# (48/p) * ( (p+8/p) - 48/p ) = 18
# (48/p) * ( p - 40/p ) = 18
# 48 - (48 * 40 / p**2) = 18
# 30 = 1920 / p**2
# p**2 = 1920 / 30 = 64 => p = 8.
# So, AH = p = 8.
# u = 48/p = 48/8 = 6. So, AG = u = 6.
# x = p + 8/p = 8 + 8/8 = 9. So, AB = AC = x = 9.
# Check with u: x = u + 18/u = 6 + 18/6 = 6 + 3 = 9. Consistent.
# So, AB = 9, AC = 9, AG = 6, AH = 8.
# GB = 9 - 6 = 3. HC = 9 - 8 = 1.
# AG*GB = 6*3=18. Correct.
# AH*HC = 8*1=8. Correct.
print("AB = AC = {}".format(9))
print("AG = {}, GB = {}".format(6, 3))
print("AH = {}, HC = {}".format(8, 1))

# Step 8: Find cos(A)
# In triangle AGH, by Law of Cosines:
# GH^2 = AG^2 + AH^2 - 2*AG*AH*cos(A)
# 5^2 = 6^2 + 8^2 - 2*6*8*cos(A)
# 25 = 36 + 64 - 96*cos(A)
# 25 = 100 - 96*cos(A)
# 96*cos(A) = 75
cos_A = 75/96
# Simplify cos_A
cos_A_num = 75
cos_A_den = 96
common_divisor = math.gcd(cos_A_num, cos_A_den)
cos_A_num //= common_divisor
cos_A_den //= common_divisor
cos_A = cos_A_num / cos_A_den
print("cos(A) = {}/{}".format(cos_A_num, cos_A_den))

# Step 9-11: Calculate AE
# Using Stewart's Theorem on triangle ABC with cevian AE:
# AC^2 * BE + AB^2 * CE = BC * (AE^2 + BE * CE)
# Since AB=AC=x, we have x^2 * (BE+CE) = BC * (AE^2 + BE*CE)
# x^2 * BC = BC * (AE^2 + BE*CE)
# x^2 = AE^2 + BE*CE
# We need BE and CE. From Menelaus: BE/CE = GB/AG = 3/6 = 1/2.
# So CE = 2*BE. BC = BE+CE = 3*BE.
# Find BC using Law of Cosines on triangle ABC:
# BC^2 = AB^2 + AC^2 - 2*AB*AC*cos(A)
# BC_squared = 9**2 + 9**2 - 2*9*9*cos_A
# BC_squared = 81 + 81 - 162 * (25/32)
# BC_squared = 162 * (1 - 25/32) = 162 * (7/32)
# BC_squared = (81*2*7)/(16*2) = (81*7)/16
# BC = sqrt(BC_squared) = 9*sqrt(7)/4
# BE = BC/3 = (3*sqrt(7))/4
# CE = 2*BE = (6*sqrt(7))/4 = (3*sqrt(7))/2
# BE_times_CE = BE * CE = (3*math.sqrt(7)/4) * (3*math.sqrt(7)/2) = (9*7)/8 = 63/8
# AE^2 = x^2 - BE*CE
# AE_squared = 9**2 - (63/8) = 81 - 63/8
# AE_squared = (648 - 63)/8 = 585/8
print("AE^2 = {}/{}".format(585, 8))

# Step 12: Express AE in the required format and find a+b+c.
# AE = sqrt(585/8) = sqrt(585) / sqrt(8) = sqrt(9*65) / (2*sqrt(2))
# AE = 3*sqrt(65) / (2*sqrt(2)) = (3*sqrt(65)*sqrt(2)) / (2*sqrt(2)*sqrt(2))
# AE = (3*sqrt(130)) / 4
# The form is a*sqrt(b)/c
a = 3
b = 130
c = 4
# Check if fraction is simplified. gcd(a,c)=gcd(3,4)=1.
# Check if b has square factors. 130 = 2*5*13. No square factors.
# The form is correct.
print("AE = ({}*sqrt({}))/{}".format(a, b, c))

# The final equation is AE = (3*sqrt(130))/4
# We need to output the numbers in the final equation.
print("So the equation is AE = (3 * sqrt(130)) / 4.")

result = a + b + c
print("a+b+c = {} + {} + {} = {}".format(a,b,c,result))