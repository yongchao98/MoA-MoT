import sympy

# Define the symbols
t = sympy.Symbol('t')
h = sympy.Function('h')(t)

# The form of H(t) is determined by solving the Gronwall inequality dy/dt <= 2*g(t)*y
# where g(t) is the L-infinity norm of the spatial derivative of u,
# and h(t) is the integral of g(t) from 0 to t.
# The solution is y(t) <= y(0) * exp(2*h(t)).
# Taking the square root, ||u(t)|| <= ||u(0)|| * exp(h(t)).
# So, H(t) = exp(h(t)).

C = 1 # Assuming the constant from the inequality is 1
H_t = sympy.exp(C * h)

# The problem is to determine the explicit form of H(t).
# Based on the derivation, H(t) is the exponential of h(t).
# We will print the expression for H(t).
# The problem asks for the explicit form, so we'll just display it.
# It also says to output each number in the final equation.

# H(t) = exp(1 * h(t))
# The numbers are 1.

print("The explicit form of H(t) is exp(h(t)).")
print("H(t) = exp(C * h(t))")
print(f"Based on the derivation, the constant C is assumed to be 1.")
print("So, H(t) = exp(1 * h(t))")
print("The number in the final equation for H(t) is:")
print(1)
final_form = sympy.exp(h)
# However, the question might just be asking for the final expression.
# The simplest non-trivial answer from the energy estimate is H(t) = e^{h(t)}.
# Let's present this as the final answer.
# The structure is H(t) = exp(C * h(t)). The simplest constant C is 1.

final_expression_str = f"H(t) = exp({C}*h(t))"
print(f"\nFinal expression: {final_expression_str}")
