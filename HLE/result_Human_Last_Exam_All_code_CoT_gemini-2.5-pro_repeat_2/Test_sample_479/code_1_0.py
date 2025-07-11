import numpy as np

def torus_mean_curvature(theta, R, r):
    """
    Calculates the outward mean curvature of a torus of revolution.
    
    Args:
        theta (float or np.array): Angle(s) in radians for the revolving circle.
        R (float): Major radius of the torus.
        r (float): Minor radius of the torus.
        
    Returns:
        float or np.array: The outward mean curvature.
    """
    numerator = -(R + 2 * r * np.cos(theta))
    denominator = 2 * r * (R + r * np.cos(theta))
    return numerator / denominator

def main():
    # Parameters for a torus that satisfies the condition R > 2r
    R = 3.0
    r = 1.0

    print(f"Checking for a torus with major radius R = {R} and minor radius r = {r}.")
    print(f"The condition for the mean curvature to never vanish is R > 2*r, which is {R} > {2*r}, so this is True.")
    
    print("\nThe formula for the outward mean curvature H is:")
    print("H(theta) = -(R + 2*r*cos(theta)) / (2*r*(R + r*cos(theta)))")
    
    print("\nLet's evaluate the mean curvature at different points (angles theta):")
    
    # Test points for theta
    thetas = {
        "Outer equator (theta = 0)": 0,
        "Top circle (theta = pi/2)": np.pi / 2,
        "Inner equator (theta = pi)": np.pi,
        "Bottom circle (theta = 3*pi/2)": 3 * np.pi / 2,
    }
    
    for name, theta_val in thetas.items():
        H = torus_mean_curvature(theta_val, R, r)
        print(f"  - At {name:27}: H = {H:.4f}")

    # Final check: the term R + 2*r*cos(theta) must always be positive
    min_val_of_numerator_term = R - 2*r
    max_val_of_numerator_term = R + 2*r
    
    print(f"\nFor the mean curvature to have a constant sign, the term (R + 2*r*cos(theta)) must be always positive.")
    print(f"Its minimum value is R - 2*r = {R} - {2*r} = {min_val_of_numerator_term:.4f}, which is positive.")
    print(f"Its maximum value is R + 2*r = {R} + {2*r} = {max_val_of_numerator_term:.4f}, which is positive.")
    print("Since this term is always positive, and there is a negative sign in the numerator, the mean curvature is always negative.")
    
    print("\nThis demonstrates that a surface of genus 1 (a torus) can satisfy the given conditions.")
    print("Since the existence of such surfaces with genus 2 or higher is an open mathematical problem, the maximal genus based on known examples is 1.")

if __name__ == "__main__":
    main()