import sympy as sp

# 1. Define the symbolic variables
# y: distance from the bottom plate
# H: height of the microchannel
# k: Debye-Huckel parameter
# z_1: base zeta potential at the bottom surface
# beta: slip length
# C1, C2: integration constants
y, H, k, z_1, beta = sp.symbols('y H k z_1 beta', real=True, positive=True)
C1, C2 = sp.symbols('C1 C2')

# 2. Define the general solution of the linearized Poisson-Boltzmann equation
# psi(y) = C1*cosh(k*y) + C2*sinh(k*y)
general_solution = C1 * sp.cosh(k * y) + C2 * sp.sinh(k * y)
psi_func = sp.Function('psi')(y)

print("The general solution is of the form:")
print(sp.Eq(psi_func, general_solution))
print("-" * 30)

# 3. Define the boundary conditions
# The slip-dependant zeta potential at the bottom wall (y=0) is z_a1
z_a1 = z_1 * (1 + beta * k)

# Boundary condition at the bottom wall: psi(0) = z_a1
bc1 = sp.Eq(general_solution.subs(y, 0), z_a1)

# Boundary condition at the top wall: psi(H) = 0
bc2 = sp.Eq(general_solution.subs(y, H), 0)

print("Applying the boundary conditions:")
print(f"1. psi(0) = z_1*(1 + beta*k)")
print(f"2. psi(H) = 0")
print("-" * 30)


# 4. Solve for the constants C1 and C2 using the boundary conditions
constants = sp.solve([bc1, bc2], (C1, C2))

# 5. Substitute the constants back into the general solution
final_solution = general_solution.subs(constants)

# 6. Simplify the resulting expression
# The sinh(A-B) identity is sinh(A)cosh(B) - cosh(A)sinh(B).
# Sympy's simplify function handles this transformation automatically.
simplified_solution = sp.simplify(final_solution)

# 7. Print the final expression for the potential distribution
print("The final expression for the Electrical Double-Layer potential distribution is:")
# The following line prints the final equation in a clear, formatted way.
# It explicitly includes the number '1' as requested by the prompt.
final_equation = sp.Eq(psi_func, simplified_solution)
print(final_equation)

# For extra clarity, let's break down the components of the final equation as text.
# The following formatting ensures that each part, including the number 1, is clearly displayed.
LHS = "psi(y)"
potential_term = f"z_1*(1 + beta*k)"
spatial_term_num = f"sinh(k*(H - y))"
spatial_term_den = f"sinh(k*H)"

print("\nIn a more descriptive format, the equation is:")
print(f"{LHS} = ( {potential_term} ) * ( {spatial_term_num} / {spatial_term_den} )")

# Final check of the equation and its components
# potential_term represents z_a1
# spatial_term represents sinh(k(H-y))/sinh(kH)
# The only explicit number in the symbolic expression is 1, which is printed.
# We will now output the final answer in the requested format.
final_answer_str = f"z_1*(1 + beta*k)*sinh(k*(H - y))/sinh(k*H)"

# The prompt asks for a single answer content, so let's format the mathematical expression.
# The user wants python code, so the printed python output is the answer.
# However, the prompt also mentions a format <<<answer content>>> at the very end.
# I will put the mathematical expression string into that format.
final_answer_formatted = f"psi(y) = {final_answer_str}"