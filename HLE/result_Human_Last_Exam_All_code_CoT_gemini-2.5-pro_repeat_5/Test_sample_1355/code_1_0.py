import math

def solve_conductance_ratio():
    """
    Calculates the ratio of the fourth moment to the average of the dimensionless
    conductance for disordered Majorana wires at the critical point.
    """

    # The problem concerns the statistical properties of conductance in disordered
    # Majorana wires at the topological critical point. This system belongs to the
    # BDI symmetry class. The dimensionless conductance 'g' follows a universal
    # distribution.
    # The n-th moment <g^n> is calculated using the formula:
    # <g^n> = [Γ(n + 1/2) * Γ(1/2)] / [2 * Γ(n + 1)]
    # where Γ is the Gamma function.

    # --- Step 1: Define a function to calculate the n-th moment ---
    def calculate_moment(n):
        # Γ(1/2) = sqrt(π)
        gamma_half = math.sqrt(math.pi)
        
        # Calculate the moment using the formula
        numerator = math.gamma(n + 0.5) * gamma_half
        denominator = 2 * math.gamma(n + 1)
        return numerator / denominator

    # --- Step 2: Calculate the average conductance <g> (n=1) ---
    # Analytically: <g> = [Γ(1.5) * Γ(0.5)] / [2 * Γ(2)]
    # Γ(1.5) = 0.5 * Γ(0.5) = 0.5 * sqrt(π)
    # Γ(2) = 1! = 1
    # <g> = [0.5 * sqrt(π) * sqrt(π)] / [2 * 1] = π / 4
    g_avg_val = calculate_moment(1)
    
    print(f"The average conductance is <g> = π / 4")
    print(f"  <g> ≈ {g_avg_val:.8f}\n")

    # --- Step 3: Calculate the fourth moment <g^4> (n=4) ---
    # Analytically: <g^4> = [Γ(4.5) * Γ(0.5)] / [2 * Γ(5)]
    # Γ(4.5) = 3.5 * 2.5 * 1.5 * 0.5 * Γ(0.5) = (105/16) * sqrt(π)
    # Γ(5) = 4! = 24
    # <g^4> = [(105/16) * sqrt(π) * sqrt(π)] / [2 * 24] = (105π / 16) / 48 = 35π / 256
    g_fourth_moment_val = calculate_moment(4)
    
    print(f"The fourth moment of the conductance is <g^4> = 35π / 256")
    print(f"  <g^4> ≈ {g_fourth_moment_val:.8f}\n")
    
    # --- Step 4: Calculate the ratio <g^4> / <g> ---
    # The ratio can be calculated from the analytical expressions, where π cancels out.
    # Ratio = (35π / 256) / (π / 4)
    num_val_ratio = 35
    den_val_ratio = 64
    final_ratio_val = num_val_ratio / den_val_ratio
    
    print("The final ratio is found by dividing the fourth moment by the average value.")
    print(f"The final equation is: <g^4> / <g> = (35π / 256) / (π / 4) = {num_val_ratio} / {den_val_ratio}")
    print(f"  Ratio = {final_ratio_val}")

solve_conductance_ratio()
<<<0.546875>>>