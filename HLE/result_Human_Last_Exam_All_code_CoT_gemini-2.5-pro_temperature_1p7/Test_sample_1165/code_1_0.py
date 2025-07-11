import sympy as sp

# Define symbols
epsilon = sp.Symbol('epsilon', positive=True, real=True)
x, s = sp.symbols('x s', real=True)
L = 1 / epsilon  # Domain length
N = L - 1      # Number of delta functions

# --- Step 1: Find the leading-order solution y0(x) ---
# The original ODE is y'' - epsilon*y' = epsilon**2 * Sum(delta(x-z_i)).
# The z_i are ordered random variables from U(0, L).
# The expected density of these points is N/L = (1/epsilon - 1) / (1/epsilon) = 1 - epsilon ≈ 1.
# The mean-field (expected) version of the ODE is y_mean'' - epsilon*y_mean' = epsilon**2.
# We propose a leading-order solution y0(x) that satisfies this mean-field equation and the boundary conditions.
# The function y0(x) = 1 - epsilon*x satisfies:
# y0(0) = 1.
# y0(L) = y0(1/epsilon) = 1 - epsilon*(1/epsilon) = 0.
# y0'' - epsilon*y0' = 0 - epsilon*(-epsilon) = epsilon**2.
# Thus, we identify y0(x) as the non-stochastic leading-order solution.
y0 = 1 - epsilon * x
print("--- Step 1: Leading-order solution y0(x) ---")
print("The leading-order (mean-field) solution satisfying the boundary conditions is:")
print("y0(x) =", y0)
print("-" * 40)

# --- Step 2: Formulate the equation for the fluctuation term Y(x) = y(x) - y0(x) ---
# Substitute y = y0 + Y into the ODE:
# (y0'' + Y'') - epsilon*(y0' + Y') = epsilon**2 * Sum(delta(x-z_i))
# Since y0'' - epsilon*y0' = epsilon**2, this simplifies to:
# Y'' - epsilon*Y' = epsilon**2 * (Sum(delta(x-z_i)) - 1)
# Y has homogeneous boundary conditions: Y(0)=0, Y(L)=0.
print("--- Step 2: Equation for the fluctuation Y(x) ---")
print("Y''(x) - epsilon*Y'(x) = epsilon**2 * (Sum(delta(x-z_i)) - 1)")
print("This is the equation for the fluctuation around the mean-field solution.")
print("-" * 40)

# --- Step 3: Analyze the stochastic forcing ---
# For small epsilon, the Y'' term dominates the left side, so: Y'' ≈ epsilon**2 * (Sum(delta(x-z_i)) - 1).
# Integrating once gives Y'(x). The core stochastic term is:
# S1(x) = Integral from 0 to x of (Sum(delta(s-z_i)) - 1) ds
# S1(x) = (Number of z_i in [0, x]) - x.
# For ordered uniform points (order statistics), S1(x) behaves statistically like a Brownian Bridge on [0, L].
# The variance of a Brownian Bridge is Var(BB(x)) = x*(L-x)/L.
# The magnitude of S1(x) fluctuations scales as sqrt(Var(S1)) ~ sqrt(L).
print("--- Step 3: Analysis of the stochastic forcing ---")
Var_S1_form = x * (L - x) / L
print("The variance of the integrated forcing term S1(x) is Var(S1(x)) ≈ x*(L-x)/L.")
print(f"The variance scales with domain size L, so the magnitude of S1(x) scales as L**(1/2) or epsilon**(-1/2).")
print("-" * 40)

# --- Step 4: Estimate the scaling of the fluctuation variance ---
# Y(x) is obtained by integrating S1(x) and applying boundary conditions.
# The variance of an integrated Brownian Bridge, Var(Integral(S1)), scales as L**3.
# So, Var(Y(x)) includes a factor of (epsilon**2)**2 from the equation and L**3 from the variance of the integrated process.
Var_Y_scaling = epsilon**4 * L**3
# Substitute L = 1/epsilon to find the final scaling with epsilon.
Var_Y_final_scaling = Var_Y_scaling.subs(L, 1/epsilon)
print("--- Step 4: Scaling of Fluctuation Variance Var(Y(x)) ---")
print(f"Var(Y(x)) ~ (epsilon**2)**2 * Var(Integral(S1)) ~ epsilon**4 * L**3")
print(f"Substituting L = 1/epsilon, Var(Y(x)) scales as: {Var_Y_final_scaling}")
print("-" * 40)

# --- Step 5: Calculate R(epsilon) ---
# The question asks for R = ( max_x|Var[y(x) - y(0)]| )^(1/2).
# Since y(0)=1 is a constant, Var[y(x) - y(0)] = Var[y(x)].
# y(x) = y0(x) + Y(x). Since y0(x) is deterministic, Var[y(x)] = Var[Y(x)].
# Therefore, R = ( max_x|Var[Y(x)]| )^(1/2).
# Since max_x Var(Y(x)) scales as epsilon, R scales as sqrt(epsilon).
R_scaling = sp.sqrt(Var_Y_final_scaling)
print("--- Step 5: Final scaling for R(epsilon) ---")
print(f"R(epsilon) ~ sqrt(max Var(Y)) ~ sqrt({Var_Y_final_scaling}) = {R_scaling}")
print("-" * 40)

# --- Step 6: Analysis of the Alternative Scenario (z_i ~ Normal(i, 0.5) i.i.d.) ---
# In this case, the z_i are independent.
# Var[Y(x)] = epsilon**4 * Var[Sum(G(x, z_i))].
# Due to independence, Var[Sum(G)] = Sum(Var[G]). For large N, this is approx. N * Var(G).
# The Green's function G(x, s) for the operator Y'' scales in magnitude with L. So Var(G) scales as L**2.
# Var(Sum(G)) ~ N * L**2.
Var_Sum_G_scaling = (N * L**2)
# Substitute N ~ L = 1/epsilon
Var_Sum_G_final = Var_Sum_G_scaling.subs(N, L)
Var_Y_iid_scaling = (epsilon**4 * Var_Sum_G_final).subs(L, 1/epsilon)
print("--- Step 6: Analysis of the i.i.d. Normal Case ---")
print(f"For independent z_i, Var[Sum(G)] scales as N * L**2 = {Var_Sum_G_final}.")
print(f"The variance Var(Y(x)) ~ epsilon**4 * (L**3), which scales as: {Var_Y_iid_scaling}.")
print("\nThe scaling for Var(Y(x)) is O(epsilon) in this case as well.")
print("Therefore, the scaling for R(epsilon) is expected to remain the same.")
print("-" * 40)