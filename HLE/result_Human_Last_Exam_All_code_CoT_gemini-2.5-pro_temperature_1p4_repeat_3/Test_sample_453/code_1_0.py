import fractions

def calculate_energy_shift():
    """
    Calculates the second-order relativistic energy shift for a hydrogen atom
    in a given state (n, l) and prints the step-by-step calculation.
    """
    # Given quantum numbers for the electron state
    n = 3
    l = 2

    print("The task is to calculate the second-order energy shift due to the relativistic kinetic energy perturbation H' = -p⁴/(8m³c²).")
    print("A direct calculation using the sum-over-states definition is prohibitively complex.")
    print("We will use the known formula for this shift derived in advanced quantum mechanics texts:")
    print("\n  ΔE = (m * c² * α⁶ / n³) * [3 / (4*n) - 2 / (2*l + 1)]\n")
    print(f"We evaluate this for the state with principal quantum number n = {n} and angular momentum quantum number l = {l}.")
    print("-" * 70)

    # --- Step 1: Calculate the term in the square brackets ---
    print("Step 1: Evaluate the term inside the square brackets: [3 / (4*n) - 2 / (2*l + 1)]")
    
    # Use the fractions module for exact rational arithmetic
    term1_denominator = 4 * n
    term2_denominator = 2 * l + 1
    term1 = fractions.Fraction(3, term1_denominator)
    term2 = fractions.Fraction(2, term2_denominator)
    
    print(f"  = [3 / (4 * {n}) - 2 / (2 * {l} + 1)]")
    print(f"  = [3 / {term1_denominator} - 2 / {term2_denominator}]")
    print(f"  = [{term1} - {term2}]")
    
    # Combine the fractions by finding a common denominator
    bracket_value = term1 - term2
    print(f"  = [({term1.numerator * term2.denominator} - {term2.numerator * term1.denominator}) / {term1.denominator * term2.denominator}]")
    print(f"  = {bracket_value}\n")

    # --- Step 2: Calculate the 1/n³ factor ---
    print(f"Step 2: Evaluate the prefactor (1/n³)")
    n_cubed_inv = fractions.Fraction(1, n**3)
    print(f"  = 1 / {n}³ = {n_cubed_inv}\n")
    
    # --- Step 3: Calculate the total numerical coefficient ---
    print("Step 3: Multiply the terms to find the total numerical coefficient.")
    total_coefficient = n_cubed_inv * bracket_value
    print(f"  Coefficient = (1 / n³) * (bracket term)")
    print(f"  = ({n_cubed_inv}) * ({bracket_value})")
    print(f"  = {total_coefficient}\n")

    # --- Step 4: Present the final symbolic answer ---
    print("-" * 70)
    print("The final result is this coefficient multiplied by the constants part of the formula (m * c² * α⁶).")
    print("\nFinal symbolic equation for the energy shift:")
    
    # Printing each number in the final equation
    print(f"ΔE = ({total_coefficient.numerator} / {abs(total_coefficient.denominator)}) * m * c² * α⁶")

if __name__ == '__main__':
    calculate_energy_shift()