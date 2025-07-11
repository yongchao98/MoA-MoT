import sympy

# Define the symbols
n = sympy.Symbol('n', positive=True, integer=True)

# Based on the derivation, the probability p_n behaves as C/n^2 for some constant C.
# p_n = C / n**2
# We want to find the limit of log(1/p_n) / log(n) as n -> infinity.

# Let's define the expression for p_n with a generic constant C.
C = sympy.Symbol('C', positive=True, constant=True)
p_n = C / n**2

# Now we compute the expression inside the limit
expression = sympy.log(1/p_n) / sympy.log(n)
simplified_expression = sympy.simplify(expression)

# The simplified expression is (log(n**2/C))/log(n) = (2*log(n) - log(C))/log(n)
# = 2 - log(C)/log(n)

# Now, we compute the limit as n tends to infinity
limit_val = sympy.limit(simplified_expression, n, sympy.oo)

# Let's print the steps for clarity
print("Let p_n be the escape probability.")
print("Based on the physical and mathematical analysis, p_n has the asymptotic form:")
print("p_n ~ C / n**2 for some constant C > 0.")
print("\nWe want to compute the limit of log(1/p_n) / log(n) as n -> infinity.")
print(f"1/p_n is proportional to n**2.")
print(f"So, log(1/p_n) is approximately log(n**2) = 2*log(n).")
print(f"The expression becomes (2*log(n)) / log(n) = 2.")
print("\nLet's formalize this with sympy:")
print(f"Expression to be evaluated: log(1/({p_n})) / log(n)")
print(f"Simplified expression: {simplified_expression}")
print(f"The limit as n -> infinity is: {limit_val}")

final_eq = f"lim_{{n->oo}} (ln(1/p_n))/(ln(n)) = lim_{{n->oo}} (ln(n^2/C))/(ln(n)) = lim_{{n->oo}} (2*ln(n) - ln(C))/(ln(n)) = 2"
print("\nThe final equation is:")
# Printing the final equation with numbers. Here, the number is 2.
print("lim (ln(1/p_n)/ln(n)) = 2")