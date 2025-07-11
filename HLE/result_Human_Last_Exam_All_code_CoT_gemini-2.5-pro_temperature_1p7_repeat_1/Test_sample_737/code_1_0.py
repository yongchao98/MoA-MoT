import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Given Data ---
    # Plate Bearing Test 1 (on subgrade)
    P_plate = 30 * 1000      # Load in N
    d_plate = 305            # Plate diameter in mm
    w1_subgrade = 2460 / 1000 # Deflection in mm

    # Plate Bearing Test 2 (on trial pavement)
    h_trial = 300            # Trial pavement thickness in mm
    w2_pavement = 1080 / 1000 # Deflection in mm

    # Design Wheel Load
    mass_design = 1.80 * 1000 # Wheel load mass in kg
    g = 9.8                  # Acceleration due to gravity in N/kg
    q_design_kpa = 600       # Tyre pressure in kN/m^2
    w_allowable = 1.00       # Max allowable deflection in mm

    # Material Properties (Poisson's Ratio mu=0.5 is assumed in the formula)
    
    print("--- Step 1: Analyze Plate Bearing Tests ---")
    
    # Calculate properties from the plate bearing test on the subgrade
    a_plate = d_plate / 2
    q_plate = P_plate / (math.pi * a_plate**2) # in N/mm^2 (MPa)
    
    # Using the formula for deflection on an elastic half-space: w = (1.5 * q * a) / E2
    # This formula is derived from the Boussinesq equation for mu=0.5
    E2_subgrade = (1.5 * q_plate * a_plate) / w1_subgrade
    print(f"Subgrade Modulus (E2) calculated from Test 1: {E2_subgrade:.2f} MPa")
    
    # For the two-layer system (Test 2), the deflection is w2 = w1 * F
    # Thus, the deflection factor F for the trial setup can be calculated directly.
    F_trial = w2_pavement / w1_subgrade
    h_a_ratio_trial = h_trial / a_plate
    print(f"Deflection Factor (F) from trial section (Test 2): {F_trial:.4f}")
    print(f"This factor corresponds to a thickness/radius ratio (h/a) of: {h_a_ratio_trial:.4f}")
    print("-" * 50)

    print("--- Step 2: Analyze Design Wheel Load ---")
    
    # Convert design pressure from kN/m^2 to N/mm^2 (MPa)
    # 1 kN/m^2 = 1000 N / (1000*1000) mm^2 = 0.001 N/mm^2
    q_design = q_design_kpa * 0.001
    
    # Calculate the total design force from the mass
    P_design = mass_design * g
    
    # Calculate the radius of the equivalent circular contact area
    a_design = math.sqrt(P_design / (math.pi * q_design))
    print(f"Design Load Force: {P_design:.0f} N")
    print(f"Design Tyre Pressure: {q_design:.3f} MPa")
    print(f"Design Contact Radius (a_design): {a_design:.2f} mm")
    print("-" * 50)
    
    print("--- Step 3: Determine Required Pavement Thickness ---")
    
    # For the design case, the deflection formula is: w_allowable = (1.5 * q_design * a_design / E2) * F_required
    # We can solve for the required deflection factor, F_required
    F_required = (w_allowable * E2_subgrade) / (1.5 * q_design * a_design)
    print(f"To limit deflection to {w_allowable:.2f} mm, the required Deflection Factor (F_required) is: {F_required:.4f}")
    
    # Core insight: Since F_required is essentially equal to F_trial, the system's
    # geometry ratio (h/a) must also be the same for the design case.
    # Therefore: h_required / a_design = h_trial / a_plate
    
    h_required = h_a_ratio_trial * a_design

    print("\nBecause F_required is the same as F_trial, we can equate the h/a ratios:")
    print("h_required / a_design = h_trial / a_plate")
    print("\nSolving for required thickness h_required:")
    print(f"h_required = (h_trial / a_plate) * a_design")
    print(f"h_required = ({h_trial:.1f} mm / {a_plate:.1f} mm) * {a_design:.1f} mm")
    print(f"The required pavement thickness is {h_required:.1f} mm.")

    return h_required

# Execute the function and print the final answer
final_thickness = solve_pavement_thickness()
# <<<190.2>>>