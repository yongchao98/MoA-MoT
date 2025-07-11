import sympy

# Define the quantum numbers
n = sympy.Integer(3)
l = sympy.Integer(2)

# Calculate intermediate values
k = l + sympy.Rational(1, 2)
l_plus_1 = l + 1

# Calculate the two terms inside the bracket of the formula
# First term: 6 * n^2 / k^2
term1_num = 6 * n**2
term1_den = k**2
term1 = term1_num / term1_den

# Second term: 8 * n^3 / (l * k * (l+1))
term2_num = 8 * n**3
term2_den = l * k * l_plus_1
term2 = term2_num / term2_den

# Calculate the pre-factor
prefactor_num = 1
prefactor_den = 16 * n**6
prefactor = sympy.Rational(prefactor_num, prefactor_den)

# Calculate the total coefficient
total_coeff = prefactor * (term1 - term2)

# Print the step-by-step calculation
print("The second-order energy shift is given by the formula:")
print("ΔE = (m*c^2*α^6 / (16*n^6)) * [ (6*n^2 / (ℓ+1/2)^2) - (8*n^3 / (ℓ*(ℓ+1/2)*(ℓ+1))) ]")
print("\nFor n = 3 and ℓ = 2:")
print(f"The pre-factor is: 1 / (16 * {n}^6) = {prefactor_num} / {prefactor_den}")
print(f"The first term in the bracket is: (6 * {n}^2) / ({l}+1/2)^2 = {term1_num} / {term1_den} = {term1}")
print(f"The second term in the bracket is: (8 * {n}^3) / ({l}*({l}+1/2)*({l}+1)) = {term2_num} / {term2_den} = {term2}")

# Print the final result
print("\nCombining these parts:")
print(f"ΔE = ({prefactor}) * ( {term1} - {term2} ) * m*c^2*α^6")
print(f"ΔE = ({prefactor}) * ( {term1 - term2} ) * m*c^2*α^6")
print(f"ΔE = {total_coeff} * m*c^2*α^6")

print("\nExpressed in terms of fundamental constants m, c, e, ħ, and ε₀, this is:")
print(f"ΔE = {total_coeff} * m * c^2 * (e^2 / (4*π*ε₀*ħ*c))^6")
