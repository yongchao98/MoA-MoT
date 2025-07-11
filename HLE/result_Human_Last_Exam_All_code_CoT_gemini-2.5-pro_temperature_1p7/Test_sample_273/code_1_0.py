import sympy as sp

# --- Step 1: Define symbolic variables ---
# lambda: scaling factor for the soliton's size
# E_A: Heisenberg exchange energy at the original size (lambda=1). Must be positive.
# E_D: Dzyaloshinskii-Moriya energy at the original size (lambda=1).
l, E_A, E_D = sp.symbols('lambda E_A E_D', real=True)
# We add the physical constraint that E_A > 0 for any non-uniform soliton.
E_A_pos = sp.Symbol('E_A', real=True, positive=True)

print("--- Scaling Analysis of 3D Soliton Stability ---")
print("The energy density is H = A*(grad(m))^2 + D*m.curl(m).")
print("We analyze the total energy E of a scaled soliton solution m(r/lambda).")
print("\nStep 1: Express the total energy E(lambda) as a function of the scaling factor lambda.")
print("In 3D, the Exchange energy (from the grad^2 term) scales as lambda^(3-2) = lambda^1.")
print("In 3D, the DMI energy (from the grad term) scales as lambda^(3-1) = lambda^2.")
# The total energy as a function of the scaling factor lambda
E_lambda = l * E_A + l**2 * E_D
print("\nTotal Energy E(lambda):")
sp.pprint(E_lambda)

# --- Step 2: Find the condition for an energy extremum (Virial Theorem) ---
print("\nStep 2: Find the condition for an energy extremum by setting the first derivative dE/dlambda = 0.")
dE_dl = sp.diff(E_lambda, l)
print("\nFirst derivative dE/dlambda:")
sp.pprint(dE_dl)

print("\nFor a static solution to exist at its natural size (lambda=1), its energy must be at an extremum:")
print("dE/dlambda | (lambda=1) = 0")
virial_eq = sp.Eq(dE_dl.subs(l, 1), 0)
print("\nThis gives the Virial Theorem for the system:")
sp.pprint(virial_eq)
print("\nWhich can be rewritten as the following equation:")
# We solve for E_A to use later
virial_solution = sp.solve(virial_eq, E_A)[0]
final_eq = sp.Eq(E_A, virial_solution)
# The prompt asks to output each number in the final equation. The numbers are 1 (implicit) and -2.
sp.pprint(final_eq)

# --- Step 3: Find the condition for stability ---
print("\nStep 3: Check the stability by ensuring the second derivative d^2E/dlambda^2 is positive.")
d2E_dl2 = sp.diff(E_lambda, l, 2)
print("\nSecond derivative d^2E/dlambda^2:")
sp.pprint(d2E_dl2)

print("\nFor the extremum to be a stable minimum, this (constant) second derivative must be positive:")
stability_condition = sp.Gt(d2E_dl2, 0)
print(f"The stability condition is: {stability_condition}")
print(f"This implies that E_D must be positive.")


# --- Step 4: Show the contradiction ---
print("\nStep 4: Combine the conditions to show the contradiction.")
print("From the Virial Theorem (Step 2): E_A = -2*E_D")
print("From the Stability Condition (Step 3): E_D > 0")

print("\nIf E_D > 0, then -2*E_D must be negative.")
print("Therefore, the Virial Theorem requires that E_A must be negative (E_A < 0).")

print("\nHowever, the exchange energy E_A = integral[A*(grad(m))^2]dV is always positive for any non-uniform soliton,")
print("since the exchange constant A and the term (grad(m))^2 are both positive.")

print("\n--- CONCLUSION ---")
print("The requirement from the stability analysis (E_A < 0) contradicts the physical nature of the exchange energy (E_A > 0).")
print("Therefore, no stable localized soliton can exist in 3D with only exchange and DMI terms.")
print("To stabilize such structures, an additional energy term with a different scaling, like a Zeeman or anisotropy term (which scale as lambda^3), is necessary.")

<<<No>>>