import math

# Step 1: Define the given constants
T = math.log(10)
# B = (1/2) * (10**20 / 99**2)
B_numerator = 10**20
B_denominator = 2 * 99**2
B = B_numerator / B_denominator

# Step 2: Use the derived relationship A^8 = 8B to find A
# This step relies on the insight that the integration limit A is the root of y_0(A)=0.
A_pow_8 = 8 * B
A = A_pow_8**(1/8)

# Step 3: Use the given relationship between A and alpha to solve for alpha
# A = (2 * alpha) / (1 - exp(-2*T))
# alpha = A * (1 - exp(-2*T)) / 2
term_in_parentheses = 1 - math.exp(-2 * T)
alpha = (A * term_in_parentheses) / 2

# Output the final result
print("The calculated value of A is:", A)
print("The calculated value of alpha is:", alpha)
print(f"""
### Final Equation Derivation ###
1. The solvability condition gives: x_0^4 + y_0^4 = K, where K = (3*alpha) / (1 - exp(-3*T)).
2. The integral is evaluated to: (K * A^4 / 4) - (A^8 / 8) = B.
3. The key simplification assumes y_0(A) = 0, which means K = A^4.
4. Substituting K=A^4 into the integral result gives: (A^8 / 4) - (A^8 / 8) = B, which simplifies to A^8 = 8*B.
5. We have A = (2*alpha) / (1 - exp(-2*T)).
6. Solving for alpha: alpha = (A * (1 - exp(-2*T))) / 2.
7. Substituting A = (8*B)^(1/8) gives the final expression for alpha.
   alpha = ((8*B)**(1/8) * (1 - exp(-2*T))) / 2
""")

T_val = math.log(10)
B_val = (1/2) * (10**20 / (99**2))

# Displaying the calculation with numbers
print("### Final Equation with Numbers ###")
print(f"alpha = ((8 * {B_val})**(1/8) * (1 - exp(-2 * {T_val}))) / 2")
print(f"alpha = (({8 * B_val})**(1/8) * ({1 - math.exp(-2*T_val)})) / 2")
print(f"alpha = ({(8 * B_val)**(1/8)} * {1 - math.exp(-2*T_val)}) / 2")
print(f"alpha = {((8 * B_val)**(1/8)) * (1 - math.exp(-2*T_val))} / 2")
print(f"alpha = {alpha}")
