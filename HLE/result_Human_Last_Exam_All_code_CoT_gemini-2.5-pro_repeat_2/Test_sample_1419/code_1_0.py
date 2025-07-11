import math

# This script explains and derives the leading order expression for the Wilson-Fisher fixed point coupling u*.

print("This script calculates the leading order expression for the Wilson-Fisher fixed point coupling u* in phi^4 theory.")
print("The theory is considered in d = 4 - epsilon dimensions, where epsilon is a small positive parameter.")

print("\n--- Step 1: The Beta Function ---")
print("The renormalization group (RG) flow describes how a coupling constant 'u' changes with the energy scale.")
print("This flow is governed by the beta function, beta(u).")
print("For the phi^4 theory with a (1/4!) * u * phi^4 interaction term, the one-loop beta function near d=4 is:")
print("beta(u) = -epsilon * u + (3 / (16 * pi^2)) * u^2 + O(u^3)")
print("We are interested in the leading order expression, so we can ignore terms of order u^3 and higher.")

print("\n--- Step 2: The Fixed Point Condition ---")
print("A fixed point, u*, is a value of the coupling where the beta function is zero.")
print("This means the coupling stops changing under the RG flow.")
print("So, we solve the equation beta(u*) = 0:")
print("-epsilon * u* + (3 * (u*)^2) / (16 * pi^2) = 0")

print("\n--- Step 3: Solving for the Fixed Point u* ---")
print("We can factor out u* from the equation: u* * (-epsilon + (3 * u*) / (16 * pi^2)) = 0")
print("This equation has two solutions:")
print("1. u* = 0, which is the trivial 'Gaussian' fixed point corresponding to a free theory.")
print("2. -epsilon + (3 * u*) / (16 * pi^2) = 0, which gives the non-trivial 'Wilson-Fisher' fixed point.")
print("\nSolving the second equation for u*:")
print("(3 * u*) / (16 * pi^2) = epsilon")
print("u* = (16 * pi^2 / 3) * epsilon")

print("\n--- Step 4: Final Expression ---")
print("The leading order expression for the fixed point coupling u* is constructed as follows.")

# Define the components of the final equation to output each number
numerator_constant = 16
denominator_constant = 3
pi_term_symbol = "pi^2"
epsilon_term_symbol = "epsilon"

print("\nThe final equation is composed of a numerator, a denominator, and symbolic terms.")
print(f"Numerator constant: {numerator_constant}")
print(f"Denominator constant: {denominator_constant}")
print(f"Pi term: {pi_term_symbol}")
print(f"Epsilon term: {epsilon_term_symbol}")

print("\nCombining these parts, we get the final expression:")
print(f"u* = ({numerator_constant} * {pi_term_symbol} / {denominator_constant}) * {epsilon_term_symbol}")

# Also, calculate the numerical value of the coefficient for reference
numerical_coefficient = numerator_constant * (math.pi**2) / denominator_constant
print(f"\nThe numerical value of the prefactor ({numerator_constant} * pi^2 / {denominator_constant}) is approximately: {numerical_coefficient:.4f}")