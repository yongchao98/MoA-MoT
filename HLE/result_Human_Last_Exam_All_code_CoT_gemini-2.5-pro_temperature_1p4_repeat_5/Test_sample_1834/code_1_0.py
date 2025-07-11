import numpy as np

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at a point P due to two
    infinite current-carrying wires.
    """
    # Define the point P where the field is to be calculated
    P = np.array([1, -1, 0])

    # For calculation purposes, we can define a constant K representing mu_0 * I / (2 * pi).
    # The actual values of mu_0 and I are not needed since they will cancel out.
    # We can set K=1.0 for this demonstration.
    K = 1.0

    # --- Calculation for Wire 1 (on x-axis, current in +x) ---

    # Perpendicular distance r1 from the x-axis to P(x,y,z) is sqrt(y^2 + z^2)
    r1 = np.sqrt(P[1]**2 + P[2]**2)
    # Magnitude of B-field from wire 1 is K / r1
    B1_mag = K / r1
    # Direction is +z (from right-hand rule)
    B1_dir = np.array([0, 0, 1])
    # B-field vector from wire 1
    B1_vec = B1_mag * B1_dir

    # --- Calculation for Wire 2 (on y-axis, current in +y) ---

    # Perpendicular distance r2 from the y-axis to P(x,y,z) is sqrt(x^2 + z^2)
    r2 = np.sqrt(P[0]**2 + P[2]**2)
    # Magnitude of B-field from wire 2 is K / r2
    B2_mag = K / r2
    # Direction is -z (from right-hand rule)
    B2_dir = np.array([0, 0, -1])
    # B-field vector from wire 2
    B2_vec = B2_mag * B2_dir

    # --- Total Magnetic Field ---

    # The total magnetic field is the vector sum (superposition)
    B_total_vec = B1_vec + B2_vec
    
    # The magnitude of the total magnetic field is the Euclidean norm of the total vector
    B_total_mag = np.linalg.norm(B_total_vec)

    # --- Print the step-by-step results ---
    
    print("This script calculates the magnetic field at P(1,-1,0).")
    print("Let K = mu_0*I/(2*pi). We can set K=1 for this calculation.")
    print("-" * 50)
    
    print("Step 1: Field from Wire 1 (on x-axis)")
    print(f"Distance r1 = sqrt({P[1]}^2 + {P[2]}^2) = {r1}")
    print(f"B1 vector = (0, 0, +K/r1) = (0, 0, {K}/{r1}) = {B1_vec}")
    print("-" * 50)

    print("Step 2: Field from Wire 2 (on y-axis)")
    print(f"Distance r2 = sqrt({P[0]}^2 + {P[2]}^2) = {r2}")
    print(f"B2 vector = (0, 0, -K/r2) = (0, 0, -{K}/{r2}) = {B2_vec}")
    print("-" * 50)

    print("Step 3: Total Magnetic Field")
    print(f"Total B vector = B1 + B2 = {B1_vec} + {B2_vec} = {B_total_vec}")
    print("-" * 50)
    
    print("Step 4: Magnitude of Total Magnetic Field")
    print("The final equation for the magnitude is: sqrt(B_total_x^2 + B_total_y^2 + B_total_z^2)")
    print(f"Magnitude = sqrt({B_total_vec[0]**2} + {B_total_vec[1]**2} + {B_total_vec[2]**2})")
    print(f"The final calculated magnitude is: {B_total_mag}")

solve_magnetic_field()