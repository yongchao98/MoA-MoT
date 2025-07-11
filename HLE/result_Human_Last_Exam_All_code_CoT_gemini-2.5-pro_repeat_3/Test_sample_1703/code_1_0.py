import math

def get_area_function_description():
    """
    This script derives and explains the formula for the area of the triangle T(t).
    It prints the components of the final equation for the area.
    """
    # --- Step 1: Define constants and derived parameters ---
    # Radius of the circle and side length of the hexagon
    s = 10.0
    # Speed of the vertices of the triangle
    v = 1.0
    
    # At t=0, the triangle T(0) is equilateral. Its vertices are midpoints of
    # alternating sides of the hexagon. The side length of this triangle, L,
    # can be shown to be 1.5 times the hexagon's side length.
    L_t0 = 1.5 * s
    
    # The area as a function of displacement 'd' from the midpoint is:
    # Area(d) = (sqrt(3)/4) * (L_t0^2 + 3 * d^2)
    # We define the coefficients for the final equation based on this.
    # Area(t) = (C1 + C2 * delta(t)^2) * sqrt(3)
    c1 = L_t0**2 / 4.0
    c2 = 3.0 / 4.0
    
    # The displacement delta(t) is periodic.
    # The vertex moves from the midpoint to one end (distance s/2), then to the other end (distance s),
    # then back to the midpoint (distance s/2).
    # Total distance for one cycle is s/2 + s + s/2 = 2*s.
    period = 2 * s / v
    
    # --- Step 2: Print the final functional form ---
    
    print("The area of the triangle T(t) is independent of the hexagon's rotation.")
    print("The area is a function of time, determined by the positions of the triangle's vertices on the hexagon's sides.")
    print("\nThe formula for the area is:\n")
    print("Area(t) = (C1 + C2 * [delta(t)]^2) * sqrt(3)\n")
    
    print("Where the numbers in the equation are:")
    print(f"  C1 = {int(L_t0**2)} / 4 = {c1}")
    print(f"  C2 = 3 / 4 = {c2}")
    print(f"  sqrt(3) is the square root of 3 (approx. {math.sqrt(3):.4f})\n")
    
    print("And delta(t) is a periodic triangular wave function representing the displacement of each vertex from the midpoint of its side:")
    print(f"  - Period: {int(period)} seconds")
    print(f"  - Amplitude: {s/2} units (half the side length of the hexagon)\n")
    
    print("The definition of delta(t) for one period (let tau = t mod 20) is:")
    print("         | tau,                 if 0 <= tau <= 5")
    print("delta(t) = | 10 - tau,            if 5 < tau <= 15")
    print("         | tau - 20,            if 15 < tau <= 20")

if __name__ == '__main__':
    get_area_function_description()