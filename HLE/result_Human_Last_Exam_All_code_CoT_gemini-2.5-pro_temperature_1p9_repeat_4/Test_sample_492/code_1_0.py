import sympy

# This script calculates the critical time 'c' for the emergence of a
# giant component in a variant Erdos-Renyi random graph model.

# Let's outline the model and the derivation step-by-step.
#
# Model:
# 1. Vertices appear at a stochastic rate of n.
# 2. Edges between any two existing vertices appear at a rate of 1/n.
# We are working in the limit where n -> infinity.
#
# Plan:
# 1. Determine the number of vertices at time t.
#    - Since vertices arrive via a Poisson process with rate n, for large n,
#      the number of vertices V(t) is approximately n*t.
#
# 2. Calculate the average degree λ(t) of the graph at time t.
#    - Consider a single vertex 'v' that arrived at time x (where 0 <= x <= t).
#    - At a later time s (where x <= s <= t), the number of other vertices is approximately n*s.
#    - The rate at which edges form connecting 'v' to other vertices at time s is
#      (number of other vertices) * (edge rate) = (n*s) * (1/n) = s.
#    - The expected degree of 'v' at time t is the integral of this rate from its
#      arrival time x to t: Integral(s ds) from x to t.
#    - The arrival times of vertices are approximately uniformly distributed in [0, t].
#      To find the average degree λ(t) of the graph, we average the expected
#      degree over all possible arrival times x in [0, t].
#      λ(t) = (1/t) * ∫[from 0 to t]( ∫[from x to t](s ds) ) dx.
#
# 3. Find the critical time c for the giant component emergence.
#    - The giant component emerges when the average degree λ exceeds 1.
#    - We set λ(c) = 1 and solve for c.

# Using sympy for symbolic mathematics
t, c, x, s = sympy.symbols('t c x s')

# --- Step 2: Derivation of the average degree λ(t) ---
print("Deriving the expression for the average degree λ(t):")

# Expected degree of a vertex that arrived at time x
degree_expr = sympy.integrate(s, (s, x, t))
# Average this over all arrival times x in [0, t]
lambda_t_expr = sympy.integrate(degree_expr, (x, 0, t)) / t

print(f"The average degree λ(t) is calculated as (1/t) * ∫[0,t] (t²/2 - x²/2) dx, which simplifies to: {lambda_t_expr}")
print("-" * 50)


# --- Step 3: Solving for the critical time c ---
print("Finding the critical time 'c' where the giant component emerges.")
print("This occurs when the average degree is 1.")

# The value the average degree must reach
val_avg_degree = 1

# Equation λ(c) = 1
equation_for_c = sympy.Eq(lambda_t_expr.subs(t, c), val_avg_degree)
print(f"We set up the equation: λ(c) = {val_avg_degree}")
print(f"The equation to solve is: c**2 / 3 = {val_avg_degree}")

# Solve the equation and show steps
num_numerator = 1
num_denominator = 3
num_rhs = 1
print(f"\nSolving the equation: c**{2} / {num_denominator} = {num_rhs}")
print(f"1. Multiply both sides by {num_denominator}:")
print(f"   c**2 = {num_rhs * num_denominator}")

print(f"2. Take the square root of both sides (we take the positive root for time):")
solutions = sympy.solve(equation_for_c, c)
positive_solution = [sol for sol in solutions if sol.is_positive][0]
print(f"   c = sqrt({num_rhs * num_denominator})")

print("-" * 50)
print("The exact value of c is the positive solution:")
# Sympy prints sqrt(3) which is the exact form.
print(f"c = {positive_solution}")

# The final answer is required in a specific format.
final_answer_str = str(positive_solution)
print(f"\n<<< {final_answer_str} >>>")