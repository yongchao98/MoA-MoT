import sympy

# Define the symbols for the calculation
# m_sq is the squared mass parameter given in the problem
# M0_sq is the squared mass of the 6th degree of freedom (the scalar mode)
# M2_sq is the squared mass of the 5 tensor modes
# a is the unknown normalization constant of the kinetic term
m_sq, M0_sq, M2_sq, a = sympy.symbols('m^2 M_0^2 M_2^2 a')

# From the problem statement, we identify the parameters:
# The mass term gives the squared mass parameter m_p^2, which we call m_sq
# The mass of the 5 DoF gives the physical squared mass M2_sq
param_eq = sympy.Eq(m_sq, m_sq)
phys_mass_eq = sympy.Eq(M2_sq, m_sq)

print("Step 1: Relate physical masses to Lagrangian parameters.")
print("The squared mass of the 5 tensor modes (M2_sq) and the 1 scalar mode (M0_sq) are related to the Lagrangian mass parameter (m_sq) and normalization (a) as follows:")
relation_M2 = sympy.Eq(M2_sq, 2 * m_sq / a)
relation_M0 = sympy.Eq(M0_sq, m_sq / a)
print(f"   {relation_M2}")
print(f"   {relation_M0}\n")

print("Step 2: Use the information from the problem statement.")
print(f"The mass of the 5 modes is given as m^2. So, we set M2_sq = m^2.")
print(f"   {phys_mass_eq}\n")

print("Step 3: Solve for the normalization constant 'a'.")
# Substitute M2_sq = m^2 into the relation for M2_sq
eq_for_a = relation_M2.subs(M2_sq, m_sq)
print(f"Substituting M2_sq with m^2 gives: {eq_for_a}")
# Solve for 'a'
a_val = sympy.solve(eq_for_a, a)[0]
print(f"Solving for 'a' gives: a = {a_val}\n")

print("Step 4: Calculate the squared mass of the sixth degree of freedom (M0_sq).")
# Substitute the value of 'a' into the relation for M0_sq
final_eq = relation_M0.subs(a, a_val)
print(f"Substituting a = {a_val} into the equation for M0_sq gives:")

# Get the components for the final print statement
# The result is M0_sq = m_sq / 2
numerator = 1
denominator = a_val
print(f"   M_0^2 = ({numerator} / {denominator}) * m^2")

# Final answer in the required format
final_answer_value = final_eq.rhs
print(f"\nThe squared mass of the sixth degree of freedom is {final_answer_value}.")