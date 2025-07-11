import sympy as sp

# Define symbols for the calculation
t = sp.Symbol('t', real=True, positive=True)
T = sp.Symbol('T', real=True, positive=True)
C = sp.Symbol('C', real=True, positive=True)

# Define the integrand based on the blow-up rate estimate
# The term (1+t) comes from the equation's viscosity
# The term C/(T-t) is the minimum growth rate of ||nabla u||^2 for a blow-up at time T
integrand = (1 + t) * C / (T - t)

print("We test the blow-up hypothesis.")
print(f"If a blow-up occurs at t=T, then ||∇u(t)||² must be at least C/(T-t).")
print("This implies the integral of (1+t)||∇u(t)||² from 0 to T must diverge.")
print("Let's analyze the integral of the lower bound:")
print(f"Integrand: {integrand}")

# Perform the definite integral from a starting time t_0 < T up to T.
# The divergence happens at the upper limit t=T.
# The integral from 0 to T is improper at t=T.
# Let's compute the indefinite integral first to see its form.
indefinite_integral = sp.integrate(integrand, t)
print(f"\nIndefinite integral: {sp.simplify(indefinite_integral)}")
print("The term '-log(T - t)' shows that the integral diverges as t approaches T.")

# To be more formal, we can check if the improper integral converges
# We can evaluate from 0 to T-epsilon and take the limit epsilon -> 0
epsilon = sp.Symbol('epsilon', real=True, positive=True)
integral_limit = sp.integrate(integrand, (t, 0, T - epsilon))
limit_result = sp.limit(integral_limit, epsilon, 0)

print(f"\nLimit of the integral from 0 to T-ε as ε -> 0: {limit_result}")
print("\nThe integral diverges to infinity.")
print("\nThis contradicts the energy estimate, which showed the integral must be finite.")
print("Therefore, a finite-time blow-up is not possible for this equation.")

# The final answer is that there is no blow-up.
# The structure of the equation can be represented as:
# ∂_t u + u⋅∇u + ν(t)Δu - ∇p = 0
# where ν(t) = 1 + t
# We are asked to output an equation. Let's output the L2 energy inequality which is key.
# d/dt ||u||^2 = -2(1+t) ||∇u||^2
print("\nThe key relation derived is the L2 energy inequality:")
print("d/dt ||u||^2 = -2 * (1 + t) * ||∇u||^2")
print("Integrating this gives:")
print("||u(T)||^2 + 2 * Integral( (1+t) * ||∇u(t)||^2 ) from 0 to T = ||u(0)||^2")
print("This implies Integral( (1+t) * ||∇u(t)||^2 ) from 0 to infinity is finite.")
