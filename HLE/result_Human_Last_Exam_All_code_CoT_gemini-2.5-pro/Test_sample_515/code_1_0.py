import math

# The problem asks for the probability that a 2D random walk, conditioned to avoid the origin,
# starting from x0 = (0,1), hits the set A of neighbors of y = (3600,0).
# The formula for this probability p(x) is approximately a(x) / a(y),
# where a(z) is the potential kernel of the 2D simple random walk.

# x0 is the starting point
x0_norm = 1.0

# y is the target point
y_norm = 3600.0

# Value of the potential kernel at a neighbor of the origin, a(x0)
# a(x0) = 2 / pi
a_x0 = 2 / math.pi

# Asymptotic value of the potential kernel for large distances ||y||
# a(y) ~= (2/pi) * log(||y||) + (2*gamma + log(8)) / pi
# where gamma is the Euler-Mascheroni constant.
gamma = 0.5772156649015328
a_y = (2 / math.pi) * math.log(y_norm) + (2 * gamma + math.log(8)) / math.pi

# The probability is the ratio a(x0) / a(y)
probability = a_x0 / a_y

# The formula can be simplified:
# prob = (2/pi) / [ (2/pi)*log(3600) + (2*gamma+log(8))/pi ]
#      = 2 / [ 2*log(3600) + 2*gamma + log(8) ]
numerator = 2
denominator = 2 * math.log(y_norm) + 2 * gamma + math.log(8)
final_prob = numerator / denominator

print("The starting point is x0 = (0,1)")
print("The target set A is the neighbors of y = (3600,0)")
print("The probability is approximately a(x0) / a(y)")
print(f"a(x0) = a((0,1)) = 2/pi = {a_x0:.4f}")
print(f"a(y) = a((3600,0)) ~= {a_y:.4f}")
print(f"Probability = {a_x0:.4f} / {a_y:.4f} = {probability:.4f}")
print("\nFinal equation:")
print(f"P = 2 / (2*log({int(y_norm)}) + 2*{gamma:.4f} + log(8))")
print(f"P = {numerator} / ({2 * math.log(y_norm):.4f} + {2*gamma:.4f} + {math.log(8):.4f})")
print(f"P = {numerator} / {denominator:.4f}")
print(f"P = {final_prob:.4f}")
print("\nApproximate answer with two significant digits:")
print(f"{final_prob:.2f}")
