import math

def calculate_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at point P(1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # Define symbolic variables for the output string
    mu_0 = "μ₀"
    I = "I"
    pi = "π"
    
    # Define the point of interest
    P = (1, -1, 0)
    
    print("This script calculates the magnetic field at P(x,y,z) = (1, -1, 0).\n")
    print("The total magnetic field is the vector sum of the fields from each wire (Principle of Superposition).")
    print("B_total = B_wire1 + B_wire2\n")
    print("The magnetic field from a single infinite wire is B = (μ₀ * I) / (2 * π * r), where r is the perpendicular distance.\n")

    # --- Calculation for Wire 1 (on x-axis) ---
    print("--- Wire 1: Current I in +x direction ---")
    r1 = math.sqrt(P[1]**2 + P[2]**2)
    print(f"The point P is at (1, -1, 0). The wire is on the x-axis.")
    print(f"The perpendicular distance r1 from the x-axis to P is sqrt(y² + z²) = sqrt(({P[1]})² + {P[2]}²) = {r1:.1f}.")
    print(f"The magnitude of the magnetic field from wire 1 is |B1| = ({mu_0} * {I}) / (2 * {pi} * {r1:.1f}).")
    print("By the right-hand rule (direction of I x direction to P), the current in the +x direction creates a field in the -z direction at P.")
    print(f"So, B1 = - ({mu_0} * {I}) / (2 * {pi}) k̂.\n")

    # --- Calculation for Wire 2 (on y-axis) ---
    print("--- Wire 2: Current I in +y direction ---")
    r2 = math.sqrt(P[0]**2 + P[2]**2)
    print(f"The point P is at (1, -1, 0). The wire is on the y-axis.")
    print(f"The perpendicular distance r2 from the y-axis to P is sqrt(x² + z²) = sqrt({P[0]}² + {P[2]}²) = {r2:.1f}.")
    print(f"The magnitude of the magnetic field from wire 2 is |B2| = ({mu_0} * {I}) / (2 * {pi} * {r2:.1f}).")
    print("By the right-hand rule (direction of I x direction to P), the current in the +y direction creates a field in the -z direction at P.")
    print(f"So, B2 = - ({mu_0} * {I}) / (2 * {pi}) k̂.\n")
    
    # --- Total Magnetic Field ---
    print("--- Total Magnetic Field ---")
    print("B_total = B1 + B2")
    print(f"B_total = [-({mu_0}*{I}) / (2*{pi})] k̂ + [-({mu_0}*{I}) / (2*{pi})] k̂")
    print(f"B_total = - 2 * (({mu_0}*{I}) / (2*{pi})) k̂")
    print(f"B_total = - ({mu_0}*{I} / {pi}) k̂\n")
    
    # --- Magnitude of the Total Field ---
    print("--- Magnitude of the Total Field ---")
    print("The magnitude is the length of the B_total vector.")
    final_magnitude_eq = f"|B_total| = ({mu_0} * {I}) / {pi}"
    print(f"Final equation for the magnitude: {final_magnitude_eq}")

if __name__ == "__main__":
    calculate_magnetic_field()
<<< (μ₀ * I) / π >>>