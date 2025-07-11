import sympy

# Define symbolic variables
rho = sympy.Symbol('rho')
pi_0 = sympy.Symbol('pi_0')
k = sympy.Symbol('k', integer=True, positive=True)

# --- Step-by-step derivation ---
print("Step 1: The detailed balance equations for a birth-death process are:")
print("pi_k * lambda_k = pi_{k+1} * mu_{k+1}")
print("\nWe can rearrange this to solve for pi_{k+1}:")
print("pi_{k+1} = (lambda_k / mu_{k+1}) * pi_k\n")

print("Step 2: Substitute the given rates:")
print("lambda_k = lambda / (k+1)")
print("mu_{k+1} = mu\n")
print("This gives:")
print("pi_{k+1} = (lambda / ((k+1) * mu)) * pi_k")
print("\nUsing the definition rho = lambda / mu, we get:")
print("pi_{k+1} = (rho / (k+1)) * pi_k\n")

print("Step 3: Express each pi_k in terms of pi_0 by iterating the relation:")
print("For k=0: pi_1 = (rho / 1) * pi_0")
print("For k=1: pi_2 = (rho / 2) * pi_1 = (rho / 2) * (rho * pi_0) = (rho**2 / 2) * pi_0")
print("For k=2: pi_3 = (rho / 3) * pi_2 = (rho / 3) * (rho**2 / 2) * pi_0 = (rho**3 / 6) * pi_0")
print("\nBy induction, the general formula is:")
print("pi_k = (rho**k / k!) * pi_0\n")

print("Step 4: Use the normalization condition that all probabilities must sum to 1:")
print("Sum_{k=0 to infinity} pi_k = 1\n")
print("Substituting our formula for pi_k:")
print("Sum_{k=0 to infinity} (rho**k / k!) * pi_0 = 1\n")
print("Factoring out pi_0:")
print("pi_0 * (Sum_{k=0 to infinity} rho**k / k!) = 1\n")

print("Step 5: Recognize the sum as the Taylor series for e^rho:")
print("The sum is the definition of e**rho:")
print("Sum_{k=0 to infinity} rho**k / k! = e**rho\n")

print("Step 6: Solve for pi_0:")
print("pi_0 * e**rho = 1")
print("Therefore, the final equation for pi_0 is:")
# Outputting the final equation as text
final_equation = sympy.Eq(pi_0, sympy.exp(-rho))
print(f"{final_equation.lhs} = {final_equation.rhs}")
