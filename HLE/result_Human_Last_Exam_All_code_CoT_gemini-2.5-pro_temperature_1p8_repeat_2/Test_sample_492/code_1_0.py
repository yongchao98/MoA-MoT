import math

# --- Introduction ---
print("This script calculates the critical time 'c' for the emergence of a giant component in a dynamic random graph.")
print("The condition for emergence is that the average vertex degree equals 1.")
print("-" * 50)

# --- Step 1: Define Model Parameters at a given time t ---
print("At any given time t, for large n:")
print(" - The number of vertices k(t) is approximately n*t.")
print(" - The rate of edge formation between any existing pair of vertices is 1/n.")
print("-" * 50)

# --- Step 2: Formulate the Average Degree d(t) ---
print("The average degree d(t) is k(t) multiplied by the average edge probability p(t).")
print("d(t) ≈ (n*t) * p(t)")
print("\nThe edge probability p(t) for two vertices depends on their arrival times T_u and T_v.")
print("p(t) ≈ E[ (1/n) * (t - max(T_u, T_v)) ]")
print("Factoring out terms gives: p(t) = (t/n) * (1 - E[max(T_u, T_v)] / t)")
print("-" * 50)

# --- Step 3: Calculate the Expected Value Term ---
# The arrival times T_u and T_v are i.i.d. Uniform(0, t).
# E[max(T_u, T_v)] / t is equivalent to E[max(X, Y)] for X, Y ~ i.i.d. Uniform(0, 1).
# This is a standard result. E[max(X, Y)] = 2/3.
E_max_XY = 2/3
print(f"The term E[max(T_u, T_v)] / t can be calculated as the expected value of the maximum of two standard uniform variables, which is {E_max_XY:.4f}.")
print("Let's derive it: if X, Y ~ U(0,1), P(max(X,Y) < z) = z^2. The PDF is 2z.")
print("E[max(X,Y)] = integral(z * 2z, z=0..1) = [2*z^3/3]_0^1 = 2/3.")
print("-" * 50)

# --- Step 4: Finalize the Expression for d(t) ---
print("Substituting this expectation back into the formulas:")
# p_t = (t/n) * (1 - E_max_XY)
# d_t = (n*t) * p_t = (n*t) * (t/n) * (1 - E_max_XY) = t^2 * (1 - 2/3) = t^2 / 3
print(f"p(t) = (t/n) * (1 - {E_max_XY:.4f}) = (t/n) * ({1-E_max_XY:.4f})")
print(f"d(t) = (n*t) * p(t) = t^2 * ({1-E_max_XY:.4f})")
print("Thus, the average degree is: d(t) = t^2 / 3")
print("-" * 50)

# --- Step 5: Solve for the Critical Time c ---
print("To find the critical time c, we set the average degree to 1:")
numerator_var = "c^2"
denominator = 3
result = 1
print(f"The final equation is: {numerator_var} / {denominator} = {result}")

# We solve for c: c^2 = 3
final_rhs = 3
print(f"Rearranging gives: c^2 = {final_rhs}")
c = math.sqrt(final_rhs)

# The exact value is sqrt(3)
exact_c = "sqrt(3)"
print(f"The exact value of c is {exact_c}")
print(f"The approximate numerical value is {c:.8f}")
