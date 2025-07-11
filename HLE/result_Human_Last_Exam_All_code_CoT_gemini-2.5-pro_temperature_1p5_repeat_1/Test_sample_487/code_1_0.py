import math

def calculate_h11(chi, c, n_fix):
    """Calculates h^{1,1}(M) based on the derived formula."""
    # The formula is h^{1,1}(M) = 13 - chi/2 + c * n_fix
    # Ensure chi is even for the division
    if chi % 2 != 0:
        return -1, "Invalid chi (must be even)"
    h11 = 13 - (chi / 2) + (c * n_fix)
    return int(h11), f"13 - ({chi} / 2) + ({c} * {n_fix}) = {int(h11)}"

def find_max_h11():
    """
    Finds the maximal h^{1,1} by checking all possible configurations of the involutions.
    """
    max_h11 = -1
    best_config_details = ""
    best_config_equation = ""

    # Set N_fix to its maximum possible value for a genus 2 curve, which is 6.
    n_fix = 6

    # --- Case 1: S^rho is a single smooth curve C_g of genus g ---
    # g is in [2, 10]
    for g in range(2, 11):
        chi = 2 - 2 * g
        c = 1  # A single curve
        h11, eq_str = calculate_h11(chi, c, n_fix)
        if h11 > max_h11:
            max_h11 = h11
            best_config_details = f"S^rho = single curve of genus g={g}"
            best_config_equation = eq_str

    # --- Case 2: S^rho = C_g U F_k (a curve and k points) ---
    # 2g - 2 + k = 8  => k = 10 - 2g. k > 0 implies g <= 4. g >= 0.
    for g in range(0, 5):
        k = 10 - 2 * g
        chi = (2 - 2 * g) + k
        c = 1 + k  # 1 curve + k points
        h11, eq_str = calculate_h11(chi, c, n_fix)
        if h11 > max_h11:
            max_h11 = h11
            best_config_details = f"S^rho = curve of genus g={g} and k={k} points"
            best_config_equation = eq_str
            
    # --- Case 3: S^rho = E_1 U E_2 (two disjoint elliptic curves) ---
    chi = 0 + 0  # Euler char of elliptic curve is 0
    c = 2      # Two components
    h11, eq_str = calculate_h11(chi, c, n_fix)
    if h11 > max_h11:
        max_h11 = h11
        best_config_details = "S^rho = two disjoint elliptic curves"
        best_config_equation = eq_str

    # --- Case 4: S^rho = E U R (an elliptic curve and a rational curve) ---
    chi = 0 + 2 # Euler char of elliptic is 0, of P^1 is 2
    c = 2       # Two components
    h11, eq_str = calculate_h11(chi, c, n_fix)
    if h11 > max_h11:
        max_h11 = h11
        best_config_details = "S^rho = one elliptic and one rational curve"
        best_config_equation = eq_str

    # --- Case 5: S^rho is a disjoint union of n rational curves and k points
    # 2g-2+k=8. g = p_a = 1-n => 2(1-n)-2+k=8 => -2n+k=8 => k=8+2n
    # chi = n*2+k = 2n+(8+2n)=4n+8
    # c = n+k = n+(8+2n)=3n+8
    # Let's check n=1 (g=0) -> k=10, chi=12, c=11. This is already in Case 2.
    
    print("The maximal value of h^{1,1}(M) is found for the following configuration:")
    print(f"Configuration: {best_config_details} on S, and N_fix={n_fix} points on C.")
    print("The calculation is:")
    print(f"h^{{1,1}}(M) = {best_config_equation}")


if __name__ == '__main__':
    find_max_h11()
