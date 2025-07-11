import sympy

# We are solving for the maximal possible genus, g.

# 1. From a theorem by Ros, for a compact surface with H < 0, the integral of the Gaussian curvature K is non-negative.
# Inequality: Integral(K dA) >= 0

# 2. From the Gauss-Bonnet theorem, the integral of the Gaussian curvature K is related to the genus g.
# Equation: Integral(K dA) = 4 * pi * (1 - g)

# 3. Combining these two results gives the inequality:
# 4 * pi * (1 - g) >= 0

# We can solve this for g. Let's represent the inequality using sympy.
g = sympy.Symbol('g')
pi = sympy.pi

# The inequality is 4*pi*(1-g) >= 0. Since 4*pi is a positive constant, we can simplify this to:
inequality = 1 - g >= 0

# Now, we solve for g
solution = sympy.solve(inequality, g)

# Let's print the steps of the final calculation.
print("The final step of the proof combines two theorems:")
print("1. Ros's Theorem implies: Integral(K dA) >= 0")
print("2. The Gauss-Bonnet Theorem states: Integral(K dA) = 4 * pi * (1 - g)")
print("\nCombining them gives the inequality:")
print("4 * pi * (1 - g) >= 0")
print("\nSince 4*pi is positive, this simplifies to:")
print("1 - g >= 0")
print("\nAdding g to both sides:")
print("1 >= g")
print(f"\nThis shows that the genus g must be less than or equal to 1. The result is: {solution}")
print("\nSince we know genus 0 (sphere) and genus 1 (torus) are possible, the maximal genus is 1.")
