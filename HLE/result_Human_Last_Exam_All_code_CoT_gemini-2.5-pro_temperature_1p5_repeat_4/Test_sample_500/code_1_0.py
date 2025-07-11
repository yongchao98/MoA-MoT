import sympy

# Introduction
print("This script determines the long-term behavior of the number of customers in the described queueing system.")
print("-" * 50)

# 1. Define the system parameters and mathematical expressions
lambda_rate = 3
u = sympy.Symbol('u', positive=True)
m = sympy.Symbol('m', positive=True, integer=True)

# The tail probability P(S > u) for large u is given.
# We check the integral of this function to determine if the mean service time E[S] is finite.
tail_integrand = 1/(3*u) + m/(u * sympy.ln(u))

print(f"Arrival rate lambda = {lambda_rate}")
print(f"The tail probability P(S > u) for large u is proportional to the expression: {tail_integrand}")
print("-" * 50)

# 2. Analyze the mean service time E[S]
# We use sympy to evaluate the integral of the tail probability from a large constant 'U0' to infinity.
U0 = sympy.Symbol('U0', positive=True)
# We need to assume U0 > 1 for ln(U0) to be positive.
integral_result = sympy.integrate(tail_integrand, (u, U0, sympy.oo))

print("Step 1: Check if the mean service time E[S] is finite.")
print(f"E[S] is given by the integral of P(S > u) from 0 to infinity.")
print(f"We check for convergence by integrating the tail of the distribution from a large value U0 to infinity:")
print(f"The result of Integral({tail_integrand}, (u, U0, oo)) is: {integral_result}")

is_finite = integral_result.is_finite
if is_finite is False:
    print("\nThe integral diverges, which means the mean service time E[S] is infinite.")
else:
    # This case will not happen based on our analysis
    print("\nThe integral converges, which means the mean service time E[S] is finite.")
print("-" * 50)

# 3. Conclude the value of liminf X_t
print("Step 2: Determine liminf X_t based on the value of E[S].")
print("The system is an M/G/infinity queue. The behavior is as follows:")
if is_finite is False:
    print("- If E[S] is infinite, the number of customers, X_t, grows to infinity almost surely.")
    print("- This implies that the limit inferior of X_t as t approaches infinity is also infinity.")
    final_answer = "infinity"
else:
    print("- If E[S] is finite, the system is stable and will return to the empty state (0 customers) infinitely often.")
    print("- In that case, the limit inferior of X_t would be 0.")
    final_answer = 0
print("-" * 50)

# There is no numerical equation for the final answer. The result is symbolic.
# We present the conclusion.
print(f"Conclusion: The value of liminf_{{t->inf}} X_t is {final_answer}.")