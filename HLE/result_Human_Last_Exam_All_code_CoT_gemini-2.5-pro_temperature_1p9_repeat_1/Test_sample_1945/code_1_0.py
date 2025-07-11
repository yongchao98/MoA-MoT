import math

# The limiting probability Pr(n) is given by the formula 9 / (2 * pi^2).
# Let's calculate and print this value.

pi = math.pi
probability = 9 / (2 * pi**2)

print("The condition for a (p,q) pair to be good is (p <= n/2) or (q <= n/2).")
print(f"The exact value of lim Pr(n) is 9 / (2 * pi^2).")
print("9 / (2 * {pi}^2) = {probability}".format(pi=pi, probability=probability))
