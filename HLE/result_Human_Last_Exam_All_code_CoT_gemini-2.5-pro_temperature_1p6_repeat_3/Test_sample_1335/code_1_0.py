# 1. Define the economic parameters from the problem description.
# All changes are in percentage points (e.g., 2% tax is input as 2.0)

# Factor cost shares
theta_LX = 2/3  # Labor's share of cost in sector X
theta_KX = 1 - theta_LX  # Capital's share of cost in sector X
theta_LY = 1/3  # Labor's share of cost in sector Y
theta_KY = 1 - theta_LY  # Capital's share of cost in sector Y

# Factor allocation shares
lambda_LX = 3/4  # Sector X's share of total labor
lambda_LY = 1 - lambda_LX  # Sector Y's share of total labor

# Elasticities
sigma_X = -2.0  # Elasticity of substitution in X
epsilon_XX = -2.0  # Own-price elasticity of demand for X
eta_IX = 1.0     # Income elasticity of demand for X

# The tax shock
tax_on_capital_X = 2.0  # 2% tax

# 2. Calculate percentage changes in wage (w) and price of good X (P_X).
# We use a "^" to denote percentage change, e.g., w_hat is the % change in w.
# Since Y is traded and capital is mobile, P_Y and r are fixed.
# Let P_Y be the numeraire, so P_Y_hat = 0.
# The rental rate for capital owners is fixed, so r_hat = 0.
# The tax applies to capital use in sector X, so r_X_hat = r_hat + tax = 2%.
# The rental rate in Y is unaffected, r_Y_hat = r_hat = 0.

print("Step 1: Calculating percentage change in nominal wage (w_hat)...")
# From the zero-profit condition in sector Y: P_Y_hat = theta_LY * w_hat + theta_KY * r_Y_hat
# 0 = (1/3) * w_hat + (2/3) * 0
w_hat = 0.0
print(f"The equation is: 0 = ({theta_LY:.2f}) * w_hat + ({theta_KY:.2f}) * 0")
print(f"Result: The percentage change in nominal wage (w_hat) is {w_hat:.2f}%\n")

print("Step 2: Calculating percentage change in the price of good X (Px_hat)...")
# From the zero-profit condition in sector X: Px_hat = theta_LX * w_hat + theta_KX * r_X_hat
Px_hat = theta_LX * w_hat + theta_KX * tax_on_capital_X
print(f"The equation is: Px_hat = ({theta_LX:.2f}) * {w_hat:.2f} + ({theta_KX:.2f}) * {tax_on_capital_X:.2f}")
print(f"Result: The percentage change in the price of good X (Px_hat) is {Px_hat:.4f}%\n")

# 3. Build and solve the general equilibrium system.

# 3a. Find production value shares (s_QX, s_QY)
# From w*L_X = theta_LX*P_X*Q_X and L_X = lambda_LX*L, we have w*lambda_LX*L = theta_LX*P_X*Q_X
# From w*L_Y = theta_LY*P_Y*Q_Y and L_Y = lambda_LY*L, we have w*lambda_LY*L = theta_LY*P_Y*Q_Y
# Dividing the two expressions: (theta_LX*P_X*Q_X)/(theta_LY*P_Y*Q_Y) = lambda_LX/lambda_LY
# This allows us to find the ratio of production values P_X*Q_X / P_Y*Q_Y
production_value_ratio = (lambda_LX / lambda_LY) * (theta_LY / theta_LX)
s_QX = production_value_ratio / (1 + production_value_ratio)
s_QY = 1 - s_QX

print("Step 3: Calculating national income and consumption changes...")
print(f"Production value shares are: s_QX = {s_QX:.2f}, s_QY = {s_QY:.2f}")

# 3b. Set up the system to solve for the change in nominal income (I_hat).
# The system links changes in consumption (C), production (Q), and income (I).
# We solve the system: I_hat * (1 - s_QX + s_QY*(lambda_LX/lambda_LY)) = s_QX*(1+eps_XX)*Px_hat - s_QY*(lambda_LX/lambda_LY)*(eps_XX*Px_hat + sigma_X*r_X_hat)
# Where r_X_hat is the tax.

lhs_I_hat_coeff = 1 - s_QX + s_QY * (lambda_LX / lambda_LY)
rhs_value = s_QX * (1 + epsilon_XX) * Px_hat - s_QY * (lambda_LX/lambda_LY) * (epsilon_XX * Px_hat + sigma_X * tax_on_capital_X)
print(f"Solving for nominal income change (I_hat) using the equation: {lhs_I_hat_coeff:.2f} * I_hat = {rhs_value:.2f}")

I_hat = rhs_value / lhs_I_hat_coeff
print(f"Result: The percentage change in nominal income (I_hat) is {I_hat:.4f}%\n")

# 3c. Calculate the change in consumption of good X (Cx_hat).
# From the demand function: Cx_hat = epsilon_XX * Px_hat + eta_IX * I_hat
print("Step 4: Calculating the final percentage change in consumption of good X (Cx_hat)...")
Cx_hat = epsilon_XX * Px_hat + eta_IX * I_hat
print(f"The equation is: Cx_hat = ({epsilon_XX:.2f}) * ({Px_hat:.4f}) + ({eta_IX:.2f}) * ({I_hat:.4f})")
print(f"Result: The percentage change in consumption of X (Cx_hat) is {Cx_hat:.4f}%\n")

# 4. Final output
print("Final calculated percentage changes:")
print("Nominal Wage:", f"{w_hat:.4f}%")
print("Price of Good X:", f"{Px_hat:.4f}%")
print("Consumption of Good X:", f"{Cx_hat:.4f}%")
print("\nFinal comma-separated values:")
print(f"{w_hat},{Px_hat},{Cx_hat}")

# For verification, converting fractions to decimals: 2/3 = 0.6666..., 29/12 = 2.4166...
# w_hat = 0.0
# Px_hat = 2/3
# Cx_hat = 29/12
final_w = 0.0
final_px = 2.0/3.0
final_cx = 29.0/12.0
# print(f"\nFinal values as fractions for precision check: 0, 2/3, 29/12")
# print(f"Which are approximately: {final_w}, {final_px}, {final_cx}")
# The values match the detailed derivation.
# The code output above will show the final result, not this check section.
# The required output format is just the numbers.

final_answer_string = f"{w_hat},{Px_hat},{Cx_hat}"
# The final answer is expected to be a string of 3 comma-separated values.
# Using the values derived: 0, 2/3, 29/12.

# Presenting the values found in the code
# w_hat = 0.0
# Px_hat = 0.666666...
# Cx_hat = 2.416666...
# Using fractions gives more precision for the final output.
final_answer_precise = f"{0.0},{2/3},{29/12}"
# The instruction is to use 'print', so the final output line will print the result string.
# To conform to the <<<...>>> format, I will output that in the last line.

# Printing the numerical results in the requested comma-separated format for the final answer.
# This must be the final line of the code output.
# print(f"\n<<< {w_hat},{Px_hat},{Cx_hat} >>>") is what i would do but I need to format it in the main text

# The main block should already print the comma separated values.
# print(f"{w_hat},{Px_hat},{Cx_hat}")

# Format the final answer string with a reasonable number of digits
# for readability while maintaining accuracy.
final_output_string = f"{w_hat:.4f},{Px_hat:.4f},{Cx_hat:.4f}"
# Let's provide this formatted string to the user.
# The user wants "0,0.6666666666666667,2.4166666666666665" or similar
# Let's not format the string to avoid losing precision
# print(final_answer_string)