import math

# The problem is to find the symbolic expression for pi_0 in terms of rho = lambda/mu.
# The following code prints the derivation steps and the final formula.

print("Step-by-step derivation of the steady-state probability pi_0:")
print("------------------------------------------------------------")

# Step 1: Define the balance equations for the birth-death process.
# The rates are lambda_i = lambda / (i + 1) and mu_i = mu.
# The balance equation is: pi_i * lambda_i = pi_{i+1} * mu_{i+1}
# Substituting the rates: pi_i * (lambda / (i + 1)) = pi_{i+1} * mu
print("1. The balance equation is pi_i * lambda_i = pi_{i+1} * mu_{i+1}.")
print("   Substituting the given rates, we get: pi_i * (lambda / (i + 1)) = pi_{i+1} * mu")

# Step 2: Express pi_n in terms of pi_0.
# From the balance equation: pi_{i+1} = pi_i * (lambda / (mu * (i + 1)))
# Let rho = lambda / mu. Then: pi_{i+1} = pi_i * (rho / (i + 1))
#
# Let's find the first few terms:
# For i=0: pi_1 = pi_0 * (rho / 1)
# For i=1: pi_2 = pi_1 * (rho / 2) = pi_0 * (rho^2 / 2!)
# For i=2: pi_3 = pi_2 * (rho / 3) = pi_0 * (rho^3 / 3!)
# The general formula is: pi_n = pi_0 * (rho^n / n!)
print("\n2. By repeatedly applying the balance equation, we find the general form for pi_n:")
print("   pi_n = pi_0 * (rho^n / n!), where rho = lambda / mu.")

# Step 3: Use the normalization condition sum(pi_n for n=0 to inf) = 1.
# sum_{n=0 to inf} [ pi_0 * (rho^n / n!) ] = 1
# pi_0 * sum_{n=0 to inf} [ rho^n / n! ] = 1
print("\n3. Using the normalization condition that all probabilities must sum to 1:")
print("   Sum(pi_n) = pi_0 * Sum(rho^n / n!) = 1")

# Step 4: Recognize the sum as the Taylor series for e^rho.
# sum_{n=0 to inf} [ rho^n / n! ] = e^rho
print("\n4. The sum is the Taylor series expansion for the exponential function e^rho.")
print("   So, the equation becomes: pi_0 * e^rho = 1")

# Step 5: Solve for pi_0.
# pi_0 = 1 / e^rho = e^(-rho)
print("\n5. Solving for pi_0 gives the final expression.")

print("\n------------------------------------------------------------")
print("Final Expression:")
# We use string formatting to display the symbolic result clearly.
# The instruction "output each number in the final equation" is interpreted
# as showing the final formula clearly.
final_equation = "pi_0 = e^(-rho)"
print(final_equation)
print("------------------------------------------------------------")

# As an example, we can calculate the value of pi_0 for a specific rho.
rho_example = 2.0
pi_0_example = math.exp(-rho_example)
print(f"\nFor instance, if rho = {rho_example}:")
print(f"pi_0 = e^(-{rho_example}) = {pi_0_example}")
