import numpy as np

def green_disk(z, zeta, R):
    """
    Computes the Green's function for the disk of radius R.
    G(z, zeta) = (1/2pi) * ln |(R - z*conj(zeta)/R) / (z - zeta)|
    """
    if np.abs(z - zeta) < 1e-9:
        raise ValueError("z and zeta must be different")

    if np.abs(zeta) < 1e-9: # Pole at the origin
        return (1/(2*np.pi)) * np.log(R / np.abs(z))
    
    # General case
    numerator = R - z * np.conj(zeta) / R
    denominator = z - zeta
    return (1/(2*np.pi)) * np.log(np.abs(numerator / denominator))

def green_regularized(p, rho, R):
    """
    Computes the regularized Green's function G(p,p) on a disk of radius rho.
    """
    if np.abs(p) < 1e-9:
        return (1/(2*np.pi)) * np.log(R/rho)
    
    return (1/(2*np.pi)) * np.log((R**2 - np.abs(p)**2) / (R * rho))

def solve():
    """
    Solves the random walk problem using the Green's function approach.
    """
    R = 1000
    z0 = 300j  # Start point (0, 300)
    a1 = 0 + 0j  # Target point (0,0)
    a2 = 2 + 0j  # Target point (2,0)
    rho = 1.0    # Effective radius for lattice points, equal to lattice spacing

    # Set up the system of linear equations to find coefficients c1, c2
    # M * c = v
    # where M is the matrix of Green's function values and v is [1, 1]
    
    M11 = green_regularized(a1, rho, R)
    M22 = green_regularized(a2, rho, R)
    M12 = green_disk(a1, a2, R)
    M21 = green_disk(a2, a1, R)
    
    M = np.array([[M11, M12], [M21, M22]])
    v = np.array([1, 1])
    
    # Solve for coefficients c = [c1, c2]
    try:
        c = np.linalg.solve(M, v)
        c1, c2 = c[0], c[1]
    except np.linalg.LinAlgError:
        print("Singular matrix, cannot solve for coefficients.")
        return

    # Calculate the probability at the starting point z0
    prob_z0 = c1 * green_disk(z0, a1, R) + c2 * green_disk(z0, a2, R)
    
    # The final equation string is for explanation. The calculation is done above.
    # The string shows the structure of the calculation with numerical values.
    # For a readable equation, let's round the intermediate values
    g_z0_a1 = green_disk(z0, a1, R)
    g_z0_a2 = green_disk(z0, a2, R)
    
    final_eq = (f"Probability = c1 * G(z0, a1) + c2 * G(z0, a2)\n"
                f"            = {c1:.3f} * {g_z0_a1:.3f} + {c2:.3f} * {g_z0_a2:.3f}\n"
                f"            = {c1 * g_z0_a1:.3f} + {c2 * g_z0_a2:.3f}\n"
                f"            = {prob_z0:.3f}")
                
    print(f"The probability that the random walk visits the set {{(0,0),(2,0)}} is calculated as follows:")
    print(final_eq)
    print("\nFinal Answer (3 significant digits):")
    print(f"{prob_z0:.3g}")
    
solve()