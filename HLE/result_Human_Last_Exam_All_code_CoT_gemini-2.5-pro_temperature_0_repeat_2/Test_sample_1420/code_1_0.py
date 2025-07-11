import sys

# Set encoding to utf-8 to ensure proper display of Greek letters
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent α using the first-order epsilon expansion.
    """
    # --- Parameters ---
    # For the Ising model universality class, the order parameter has n=1 component.
    n = 1
    # The spatial dimension of the system.
    d = 3
    # The upper critical dimension for this universality class.
    d_c = 4

    print("Calculating the specific heat scaling exponent α for the Ising model (n=1) in d=3.")
    print("The calculation uses the epsilon (ϵ) expansion to first order.")
    print("-" * 60)

    # --- Step 1: Calculate epsilon (ϵ) ---
    # Epsilon is the deviation from the upper critical dimension d_c.
    epsilon = d_c - d
    print(f"Step 1: Calculate the expansion parameter ϵ = d_c - d")
    print(f"ϵ = {d_c} - {d} = {epsilon}")
    print("")

    # --- Step 2: Apply the formula for α ---
    # The first-order epsilon expansion formula for α is:
    # α = (4 - n) / (2 * (n + 8)) * ϵ
    print("Step 2: Apply the first-order epsilon expansion formula for α.")
    print("Formula: α = (4 - n) / (2 * (n + 8)) * ϵ")
    print(f"Substituting n = {n} and ϵ = {epsilon} into the formula:")

    # Perform the calculation step-by-step
    numerator = 4 - n
    n_plus_8 = n + 8
    denominator = 2 * n_plus_8
    alpha = (numerator / denominator) * epsilon

    print(f"α = ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator} / (2 * {n_plus_8}) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}")
    print("")

    # --- Final Result ---
    print(f"The final calculated value for the scaling exponent α is: {alpha}")
    return alpha

if __name__ == "__main__":
    calculate_alpha_exponent()