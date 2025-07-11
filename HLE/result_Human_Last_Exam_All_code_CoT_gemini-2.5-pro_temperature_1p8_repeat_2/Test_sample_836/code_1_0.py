import sympy
from sympy import Function, Integral, Symbol, exp, I, oo, Eq, pprint
from sympy.vector import CoordSys3D, dot

def display_fourier_transform_definition():
    """
    This function symbolically defines and prints the space-time double Fourier
    transform of the generalized pair correlation function.
    """

    # --- Introduction and Explanation ---
    print("In the nuclear criticality community, the space-time double Fourier transform")
    print("of the generalized pair correlation function G(r, t) is called the")
    print("Dynamic Structure Factor, commonly denoted as S(k, ω).\n")
    
    print("1. The Generalized Pair Correlation Function, G(r, t):")
    print("   This function describes the correlation in the fluctuations of the neutron")
    print("   density. Specifically, it measures the likelihood of finding a neutron at a")
    print("   certain position and time relative to another neutron at a reference position")
    print("   and time. It is a function of the spatial separation (r) and time lag (t)")
    print("   between the two points.\n")

    print("2. The Dynamic Structure Factor, S(k, ω):")
    print("   This is the result of Fourier transforming G(r, t) with respect to both")
    print("   space (r) and time (t). The new variables are the wavevector (k) and")
    print("   angular frequency (ω). S(k, ω) reveals the spectrum of the density")
    print("   fluctuations, showing which wavelengths (related to 1/k) and frequencies (ω)")
    print("   are dominant in the system's neutron noise.\n")

    print("--- Mathematical Definition ---")
    print("The symbolic equation defining this relationship is S(k, ω) = ∫∫ G(r, t) * e^(i*(k·r - ωt)) d³r dt:")

    # --- Symbolic Representation using SymPy ---
    
    # Define time and frequency variables
    t = Symbol('t', real=True, doc="Time lag")
    omega = Symbol('ω', real=True, doc="Angular frequency")

    # Define a 3D Cartesian coordinate system for space and wavevector
    # 'r' is the position vector (x, y, z)
    # 'k' is the wavevector (k_x, k_y, k_z)
    R_sys = CoordSys3D('R')
    r_vec = R_sys.x * Symbol('x') + R_sys.y * Symbol('y') + R_sys.z * Symbol('z')
    k_vec = R_sys.x * Symbol('k_x') + R_sys.y * Symbol('k_y') + R_sys.z * Symbol('k_z')
    
    # Shorthand for volume element d³r = dx dy dz
    d3_r = (R_sys.x.name, -oo, oo), (R_sys.y.name, -oo, oo), (R_sys.z.name, -oo, oo)

    # Define the functions G and S
    G = Function('G')(R_sys.x, R_sys.y, R_sys.z, t)
    S = Function('S')(k_vec, omega)

    # Construct the integrand based on the physics convention
    # Note: exp(i*(k·r - ωt)) is a common convention
    integrand = G * exp(I * (dot(k_vec, r_vec) - omega * t))

    # Create the integral expression over all space (d³r) and all time (dt)
    transform = Integral(integrand, d3_r, (t, -oo, oo))

    # Form the complete equation S(k, ω) = Integral(...)
    final_equation = Eq(S, transform, evaluate=False)
    
    # Pretty-print the final equation
    print() # Add a newline for better formatting
    pprint(final_equation, use_unicode=True)
    print()


if __name__ == '__main__':
    display_fourier_transform_definition()