import math

def calculate_conductance_moment(n):
    """
    Calculates the n-th moment of the dimensionless conductance 'g' for a
    disordered Majorana wire at the critical point.

    The formula is: <g^n> = Gamma(n + 1/2) / (2 * sqrt(pi) * Gamma(n + 1))
    """
    numerator = math.gamma(n + 0.5)
    denominator = 2 * math.sqrt(math.pi) * math.gamma(n + 1)
    return numerator / denominator

def main():
    """
    Calculates the ratio between the fourth statistical moment of the
    dimensionless conductance and its average value.
    """
    print("This script calculates the ratio <g^4> / <g> for a disordered Majorana wire at the critical point.")
    print("The moments are given by the formula: <g^n> = Gamma(n + 1/2) / (2 * sqrt(pi) * Gamma(n + 1))\n")

    # --- Step 1: Calculate the average conductance <g> (n=1) ---
    n_avg = 1
    avg_g = calculate_conductance_moment(n_avg)
    # The exact value is 1/4
    
    print(f"Step 1: Calculate the average conductance <g> for n = {n_avg}")
    print(f"   <g> = {avg_g:.8f}  (Exact value: 1/4)")
    print("-" * 50)

    # --- Step 2: Calculate the fourth moment <g^4> (n=4) ---
    n_fourth = 4
    fourth_moment_g = calculate_conductance_moment(n_fourth)
    # The exact value is 35/256

    print(f"Step 2: Calculate the fourth moment <g^4> for n = {n_fourth}")
    print(f"   <g^4> = {fourth_moment_g:.8f} (Exact value: 35/256)")
    print("-" * 50)

    # --- Step 3: Compute the final ratio ---
    if avg_g == 0:
        print("Error: Average conductance is zero, cannot compute ratio.")
        return
        
    ratio = fourth_moment_g / avg_g
    # The exact value is (35/256) / (1/4) = 35/64
    
    print("Step 3: Compute the final ratio <g^4> / <g>")
    # As requested, showing the numbers in the final equation
    print(f"   Ratio = <g^4> / <g> = {fourth_moment_g:.8f} / {avg_g:.8f}")
    print(f"\n   The final result is: {ratio:.8f}")
    print("   The exact fractional value is 35/64.")

if __name__ == "__main__":
    main()
