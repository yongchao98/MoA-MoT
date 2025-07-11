print("The deduced system of equations for the baseline simulation (Simulation 1) is:")

# Parameters deduced from plot analysis:
# k=2, so Re = 100.
# b = -1
# c = 10
# d = -1
# e = -1
# Parameters a1, a2 remain undetermined from the set {-1, 10}.
# We will represent them symbolically.

print("x1' = -2*x1 + a1*x2*x3 + a2*x4*x5")
print("x2' = -9*x2 - 1*x1*x3")
print("x3' = -5*x3 - 10*x1*x2 + 100")
print("x4' = -5*x4 + 1*x1*x5")
print("x5' = -x5 + 1*x1*x4")

print("\nOutputting each number in the final equations as requested:")
print("Equation 1 numbers: -2 (a1 and a2 are undetermined parameters from {-1, 10})")
print("Equation 2 numbers: -9, -1")
print("Equation 3 numbers: -5, -10, 100")
print("Equation 4 numbers: -5, 1")
print("Equation 5 numbers: -1, 1")