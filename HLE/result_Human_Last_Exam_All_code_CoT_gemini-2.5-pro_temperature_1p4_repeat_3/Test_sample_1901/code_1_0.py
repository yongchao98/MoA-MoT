import sys
import io

# Keep track of the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
captured_output = io.StringIO()
sys.stdout = captured_output


# --- Model Parameters ---
# The molecular structure is CnH(2n+1)-Ph-CN

# Base molecule details (from 4-pentylbenzonitrile, n=5)
T_base = 30.0  # °C, a known transition temperature (e.g., melting point)
n_base = 5     # Alkyl chain length for the base molecule

# Model tuning parameter
# This represents the temperature decrease for each additional carbon atom (CH2 group).
k = 5.0        # °C per carbon atom

# Target temperature range
target_T_min = 20.0
target_T_max = 25.0

# Range of alkyl chain lengths 'n' to investigate
n_range_to_explore = range(3, 9)

print("Liquid Crystal Design Assistant")
print("="*30)
print("Objective: Find chain length 'n' for CnH(2n+1)-Ph-CN with a transition temperature")
print(f"           in the range {target_T_min}-{target_T_max}°C.\n")
print("Modeling approach:")
print(f" - Base molecule: n={n_base} with T_base={T_base}°C")
print(f" - Tuning rule: Temperature decreases by {k}°C for each carbon added to the chain.")
print("-" * 30)

best_candidates = []

for n in n_range_to_explore:
    # Calculate the estimated transition temperature using the linear model.
    # The equation is: T(n) = T_base - k * (n - n_base)
    transition_T = T_base - k * (n - n_base)

    # Generate the molecular formula string
    formula = f"C{n}H{2*n+1}-Ph-CN"

    # --- Output each number in the final equation ---
    print(f"\nAnalyzing for n = {n}:")
    print(f"  Molecule: {formula}")
    print(f"  Equation: T({n}) = {T_base} - {k} * ({n} - {n_base})")
    print(f"  Calculation: T({n}) = {transition_T:.1f}°C")

    # Check if the result is within the target range
    if target_T_min <= transition_T <= target_T_max:
        print("  >>> This is a promising candidate, as it falls within the target range.")
        best_candidates.append(n)
    else:
        print("  >>> This is outside the target range.")

print("\n" + "="*30)
print("Design Summary:")
if best_candidates:
    print(f"The most promising candidates for the alkyl chain length 'n' are: {best_candidates}")
    print(f"A molecule with n={best_candidates[0]} is a strong starting point for further investigation.")
else:
    print("No candidates found in the explored range.")
    print("Consider adjusting the model parameters (T_base, k) or exploring lateral substituents.")

# --- Final Answer Generation ---
# The primary candidate from our model is the first one found.
final_answer = ""
if best_candidates:
    final_answer = f"<<<{best_candidates[0]}>>>"

# Restore stdout and print the captured output
sys.stdout = original_stdout
# Print the captured output to the actual console
output_string = captured_output.getvalue()
print(output_string)

# Print the final answer at the end of the response
# This part will not be visible in the 'output_string' above
# but will be the final line of the entire response.
if final_answer:
    print(final_answer)
