import math

# The problem asks for the probability that a 2D random walk conditioned to avoid the origin,
# starting from x0 = (0,1), hits the set of neighbors of z = (3600,0).
# This probability can be approximated by the ratio of the potential kernel values a(x0)/a(z).

# The potential kernel a(x) for the 2D simple random walk has known values and asymptotics:
# 1. For a neighbor of the origin, like (0,1), the value is exact: a((0,1)) = 4/pi.
# 2. For a point x far from the origin, a(x) is approximately:
#    a(x) ~ (2/pi) * ln(||x||) + (1/pi) * (gamma + ln(8))
#    where gamma is the Euler-Mascheroni constant.

# Let x0 = (0,1) and z = (3600,0). The probability P is approximately a(x0) / a(z).
# P ~ (4/pi) / [ (2/pi)*ln(||z||) + (1/pi)*(gamma + ln(8)) ]
# P ~ 4 / [ 2*ln(3600) + gamma + ln(8) ]
# P ~ 4 / [ 2*ln(3600) + gamma + 3*ln(2) ]

# Constants
R = 3600.0
# Euler-Mascheroni constant
gamma = 0.57721566490153286060651209008240243104215933593992

# Numerator of the expression for the probability
numerator = 4

# Denominator of the expression for the probability
denominator = 2 * math.log(R) + gamma + 3 * math.log(2)

# Calculate the probability
probability = numerator / denominator

print("The probability is calculated using the formula P = a(x0) / a(z)")
print("where x0 = (0,1) and z = (3600,0).")
print("P approx= 4 / (2*ln(3600) + gamma + 3*ln(2))")
print(f"P approx= {numerator} / (2*ln({R}) + {gamma:.4f} + 3*ln(2))")
print(f"P approx= {numerator} / ({2 * math.log(R):.4f} + {gamma:.4f} + {3 * math.log(2):.4f})")
print(f"P approx= {numerator} / {denominator:.4f}")
print(f"P approx= {probability:.2g}")

# The final answer should be given with two significant digits.
# The format is <<<answer>>>
final_answer = f"{probability:.2g}"
# print(f"<<<{final_answer}>>>")