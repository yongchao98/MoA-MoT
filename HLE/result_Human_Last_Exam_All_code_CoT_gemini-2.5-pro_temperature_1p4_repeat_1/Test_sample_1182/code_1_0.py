import math

def get_poly_str(coeffs):
    """Creates a string representation of a polynomial from its coefficients."""
    poly_str = []
    for i, c in reversed(list(enumerate(coeffs))):
        if c == 0:
            continue
        if i == 0:
            poly_str.append(f"{c}")
        elif i == 1:
            poly_str.append(f"{c}*x")
        else:
            poly_str.append(f"{c}*x^{i}")
    return " + ".join(poly_str).replace("+ -", "- ")

def get_degree(coeffs):
    """Gets the degree of a polynomial."""
    for i in range(len(coeffs) - 1, -1, -1):
        if coeffs[i] != 0:
            return i
    return -1

def calculate_genus_hyperelliptic(coeffs):
    """Calculates the genus of a hyperelliptic curve y^2=f(x)."""
    deg = get_degree(coeffs)
    if deg < 0:
        return 0
    return (deg - 1) // 2

def get_genus_reg_F2(coeffs_F2):
    """
    Calculates the genus of the curve y^2 = f(x) over F_2.
    """
    deg = get_degree(coeffs_F2)
    if deg <= 2:
        # y^2 = ax^2 + bx + c -> (y + sqrt(a)x)^2 = ... is genus 0
        return 0
    # For y^2 = f(x) over F_2, if f is not a square
    # the genus is given by the standard formula.
    return (deg - 1) // 2


def solve():
    """
    Finds the number of double points in the stable reduction of the given curve.
    """
    # Coefficients of f(x) = 8*x^5 + 4*x^4 + 4*x^3 + 1*x^2 + 8*x + 0
    # from lowest degree (x^0) to highest (x^5)
    coeffs = [0, 8, 1, 4, 4, 8]
    
    print(f"Original curve: y^2 = {get_poly_str(coeffs)}")

    # Calculate genus of the original curve
    g_orig = calculate_genus_hyperelliptic(coeffs)
    
    # Loop to find the regular model
    for i in range(10): # Limit iterations to prevent infinite loops
        # Reduce coefficients modulo 2
        coeffs_F2 = [c % 2 for c in coeffs]
        
        # Check if the reduced polynomial is a square in F_2.
        # This is true if all odd-degree coefficients are 0.
        is_square = True
        for j in range(1, len(coeffs_F2), 2):
            if coeffs_F2[j] != 0:
                is_square = False
                break
        
        # If it is not a square, we have found our regular model's reduction.
        if not is_square:
            print(f"\nAfter {i} transformations, the model reduction is not a square.")
            break
            
        print(f"\nTransformation {i+1}:")
        
        # Apply the transformation: c_i -> c_i * 2^(i-2)
        new_coeffs = [0] * len(coeffs)
        for j in range(len(coeffs)):
            if j >= 2:
                new_coeffs[j] = coeffs[j] * (2**(j - 2))
            else:
                # Use integer division for negative exponents
                new_coeffs[j] = coeffs[j] // (2**(2 - j))
        coeffs = new_coeffs
        print(f"New model: y^2 = {get_poly_str(coeffs)}")
    
    # The final regular model
    print(f"\nThe equation for the final regular model is:\ny^2 = {get_poly_str(coeffs)}")
    
    # Reduction mod 2 of the final model
    final_coeffs_F2 = [c % 2 for c in coeffs]
    print(f"\nIts reduction modulo 2 is:\ny^2 = {get_poly_str(final_coeffs_F2)}")

    # Genus of the reduced curve
    g_reg = get_genus_reg_F2(final_coeffs_F2)
    
    # The number of double points is the difference in genera
    delta = g_orig - g_reg
    
    print("\n--- Calculation ---")
    print(f"Genus of the original curve: g = {g_orig}")
    print(f"Genus of the reduced curve's smooth model: g_reg = {g_reg}")
    print(f"The number of double points in the stable reduction is delta = g - g_reg")
    print(f"delta = {g_orig} - {g_reg} = {delta}")
    print("\nFinal answer:")
    print(f"{delta}")

solve()
<<<2>>>