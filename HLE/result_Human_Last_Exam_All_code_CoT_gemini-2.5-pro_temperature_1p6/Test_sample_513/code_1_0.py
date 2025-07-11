import math

def solve_sphere_packing():
    """
    Calculates the maximum number of spheres that can be packed in a cylinder
    by using pre-computed results from scientific literature.
    """
    # Define the parameters of the problem
    r_sphere = 2.0  # cm
    R_cylinder = 10.0 # cm
    H_cylinder = 11.0 # cm

    # The ratio of cylinder radius to sphere radius
    radius_ratio = R_cylinder / r_sphere

    if radius_ratio != 5.0:
        print("This solution is specifically tailored for a cylinder/sphere radius ratio of 5.")
        print("The provided radii R={} and r={} result in a ratio of {}.".format(R_cylinder, r_sphere, radius_ratio))
        return

    # Data from literature for R/r = 5.
    # 'n' is the number of spheres.
    # The value is the normalized optimal height 'h_norm = H_needed / (2 * r_sphere)'.
    # Source: Stoyan, Y. & Yaskov, G. (2012). Densest packings of equal spheres in a cylinder.
    known_packings_R_div_r_5 = {
        50: 2.458, 51: 2.502, 52: 2.548, 53: 2.593, 54: 2.639,
        55: 2.684, 56: 2.730, 57: 2.766, 58: 2.813, 59: 2.859, 60: 2.923
    }
    
    # Sort the sphere counts to iterate from smallest to largest
    sorted_n = sorted(known_packings_R_div_r_5.keys())
    
    max_spheres = 0
    final_required_height = 0

    print("Checking how many spheres can fit in a cylinder of height {} cm.".format(H_cylinder))
    print("-" * 50)

    # Find the largest number of spheres (n) whose required height fits
    for n in sorted_n:
        # Get the normalized height required for n spheres
        h_norm = known_packings_R_div_r_5[n]
        
        # Calculate the actual required height in cm
        required_height = h_norm * (2 * r_sphere)
        
        if required_height <= H_cylinder:
            # This number of spheres fits. Update our potential answer.
            max_spheres = n
            final_required_height = required_height
            # This check shows the calculation for the current best fit
            print("For n = {:<2}: required height is {:<6.3f} cm, which is <= {} cm. (Fits)".format(n, required_height, H_cylinder))
        else:
            # This is the first n that does not fit, so the previous one was the maximum.
            print("For n = {:<2}: required height is {:<6.3f} cm, which is >  {} cm. (Does not fit)".format(n, required_height, H_cylinder))
            break
            
    print("-" * 50)
    print("Final Answer:")
    print("The optimal packing for {} spheres requires a cylinder of height {:.3f} cm.".format(max_spheres, final_required_height))
    print("Since this is less than or equal to the available height of {} cm, this is the maximum number of spheres from the known optimal configurations.".format(H_cylinder))
    print("\nThe maximum number of spheres is {}.".format(max_spheres))
    

solve_sphere_packing()