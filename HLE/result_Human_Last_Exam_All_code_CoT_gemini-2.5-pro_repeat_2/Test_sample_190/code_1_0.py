import sympy

# Define symbolic variables for the state k and the parameter c
k, c = sympy.symbols('k c')

# --- Introduction and Problem Setup ---
print("We want to find the infimum of c for which the given Markov chain is transient.")
print("The transition probabilities for large state k are:")
print(f"P(k -> k-2) = 1/4")
print(f"P(k -> k+2) = 1/4")
print(f"P(k -> k-1) = 1/4 - c/k")
print(f"P(k -> k+1) = 1/4 + c/k")
print("-" * 50)

# --- Step 1: Calculate Mean Drift ---
print("Step 1: Calculate the expected drift mu_k = E[X_{n+1} - k | X_n = k].")

# Define transition probabilities using sympy for precise calculations
p_k_minus_2 = sympy.Rational(1, 4)
p_k_plus_2 = sympy.Rational(1, 4)
p_k_minus_1 = sympy.Rational(1, 4) - c/k
p_k_plus_1 = sympy.Rational(1, 4) + c/k

# Jumps from state k
jumps = [-2, -1, 1, 2]
probs = [p_k_minus_2, p_k_minus_1, p_k_plus_1, p_k_plus_2]

# Calculate the expected drift
mu_k = sum(jump * prob for jump, prob in zip(jumps, probs))
mu_k_simplified = sympy.simplify(mu_k)

print(f"mu_k = (-2) * ({p_k_minus_2}) + (-1) * ({p_k_minus_1}) + (1) * ({p_k_plus_1}) + (2) * ({p_k_plus_2})")
print(f"After simplification, mu_k = {mu_k_simplified}")

# The drift is of the form A/k. We identify the coefficient A.
A = mu_k_simplified * k
print(f"This drift is of the form A/k, where A = {A}\n")

# --- Step 2: Calculate Asymptotic Second Moment ---
print("Step 2: Calculate the asymptotic second moment B = lim_{k->inf} E[(X_{n+1} - k)^2 | X_n = k].")

# Squared jumps from state k
squared_jumps = [j**2 for j in jumps]

# Calculate the second moment of the jump
m2_k = sum(sq_jump * prob for sq_jump, prob in zip(squared_jumps, probs))
m2_k_simplified = sympy.simplify(m2_k)

print(f"M2_k = E[(X_{n+1} - k)^2] = (-2)^2*({p_k_minus_2}) + (-1)^2*({p_k_minus_1}) + (1)^2*({p_k_plus_1}) + (2)^2*({p_k_plus_2})")
print(f"After simplification, M2_k = {m2_k_simplified}")

# B is the limit of M2_k as k -> infinity
B = sympy.limit(m2_k_simplified, k, sympy.oo)
print(f"The limit as k -> infinity gives B = {B}\n")

# --- Step 3: Apply Transience Criterion ---
print("Step 3: Apply the transience criterion A > B/2.")

# The inequality for transience is A > B/2
inequality_rhs = B / 2

print("The final equation to determine transience is:")
print(f"    A > B / 2")
print("Substituting the values we found for A and B:")
# Outputting each number in the final equation as requested
print(f"    {A} > {B} / 2")
print(f"    {A} > {inequality_rhs}")

# --- Step 4: Find the Infimum ---
print("\nStep 4: Solve for c to find the infimum.")

# We solve the inequality 2*c > 5/4 for c
# c > 5/8
solution = sympy.solve(sympy.Gt(A, inequality_rhs), c)

print(f"Solving the inequality for c gives: {solution}.")
infimum_value = inequality_rhs / (A/c)
print(f"Thus, the chain is transient if c > {infimum_value}.")
print(f"The infimum of this set of c values is {infimum_value}.")

<<<5/8>>>