import numpy as np

def potential_kernel(z):
    """
    Calculates the approximate potential kernel a(z) for the 2D SRW
    using its asymptotic formula for a point z = (z1, z2).
    """
    z_vec = np.array(z)
    norm_z = np.linalg.norm(z_vec)
    if norm_z == 0:
        return 0
    # The constant is C = (2/pi) * (gamma + log(8)), where gamma is the Euler-Mascheroni constant.
    # C is approximately 1.693.
    C = (2 / np.pi) * (np.euler_gamma + np.log(8))
    return (2 / np.pi) * np.log(norm_z) + C

def green_function_h_transform(x, y):
    """
    Calculates the Green's function G_hat(x, y) for the h-transformed process.
    """
    x_vec = np.array(x)
    y_vec = np.array(y)
    
    ax = potential_kernel(x_vec)
    ay = potential_kernel(y_vec)
    axy = potential_kernel(x_vec - y_vec)
    
    if ax == 0:
        # Avoid division by zero if x is the origin, although the problem is on Z^2\{0}
        return float('inf')
        
    # G_0(x, y) = a(x) + a(y) - a(x-y)
    G0_xy = ax + ay - axy
    
    # G_hat(x, y) = (a(y) / a(x)) * G_0(x, y)
    G_hat_xy = (ay / ax) * G0_xy
    return G_hat_xy

def main():
    # Fixed starting point x
    x = (2, 0)

    print("Demonstrating the growth of the Green's function G_hat(x, y) as ||y|| increases.")
    print(f"Fixed starting point x = {x}\n")

    # Sequence of points y moving away from the origin along the x-axis
    y_points = [(k, 0) for k in range(10, 101, 10)]

    for y in y_points:
        g_hat_value = green_function_h_transform(x, y)
        print(f"For y = {str(y):<8}, G_hat({x}, {y}) = {g_hat_value:.4f}")

if __name__ == "__main__":
    main()