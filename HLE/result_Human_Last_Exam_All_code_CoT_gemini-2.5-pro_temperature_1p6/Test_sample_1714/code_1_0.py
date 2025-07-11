def solve_pandora_mystery():
    """
    This script calculates the dark matter percentage of the Pandora galaxy,
    simulating the constraints of the Wuxing computing architecture.

    The Wuxing system does not support floating point numbers, relying on a
    special `frac` type. Key constants must be approximated to fit within
    the `frac` type's limited-size numerator (signed char, -128 to 127).

    - Gravitational Constant (G) in galaxy-friendly units is ~4.302e-6.
      We approximate this to 4e-6 to be representable.
    - Other values are given: v=200 km/s, R=10 kpc, M_lum=6e9 M_sun.
    """

    # --- Step 1: Define constants and inputs based on Wuxing constraints ---

    # Approximated Gravitational Constant G in units of (kpc * (km/s)^2 / M_sun)
    G_const = 4e-6

    # Given data
    v = 200  # velocity in km/s
    R = 10   # radius in kpc
    M_lum = 6e9  # Luminous mass in Solar Masses

    # --- Step 2: Perform the calculations ---

    # Total Mass M_total = (v^2 * R) / G
    M_total = (v**2 * R) / G_const

    # Dark Mass = Total Mass - Luminous Mass
    M_dark = M_total - M_lum

    # Percentage = (Dark Mass / Total Mass) * 100
    percentage = (M_dark / M_total) * 100
    
    # --- Step 3: Output the results as requested ---
    
    # The final equation is percentage = (M_total - M_lum) / M_total * 100
    # The prompt requires printing each number in this final equation.
    print("This script simulates the calculation on the Wuxing architecture.")
    print("Values used in the final equation:")
    print(f"Total Mass (M_total) = {M_total:.1e} M_sun")
    print(f"Luminous Mass (M_lum) = {M_lum:.1e} M_sun")
    print("\nFinal Equation Calculation:")
    print(f"({M_total:.1e} - {M_lum:.1e}) / {M_total:.1e} * 100 = {percentage:.1f}%")

solve_pandora_mystery()