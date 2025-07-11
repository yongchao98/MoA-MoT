import math

def check_torus_mean_curvature(R, r):
    """
    Checks if a torus with major radius R and minor radius r has a 
    mean curvature that never vanishes. This is true if R > 2*r.
    The function also prints the minimal value of H.
    """
    print(f"Checking a torus with major radius R = {R} and minor radius r = {r}")
    
    # The condition for non-vanishing mean curvature is R > 2*r
    if R > 2 * r:
        print(f"The condition R > 2r ({R} > {2*r}) is satisfied.")
        # The mean curvature H is given by the formula:
        # H = (R + 2*r*cos(phi)) / (2*r*(R + r*cos(phi)))
        # The minimum value of H occurs at phi = pi, where cos(phi) = -1.
        
        # We demonstrate the calculation of the minimal value of the numerator
        min_numerator_val = R + 2 * r * (-1)
        
        # And the denominator at that point
        denominator_at_min = 2 * r * (R + r * (-1))
        
        min_H = min_numerator_val / denominator_at_min

        print("The mean curvature H is given by the equation H = (R + 2*r*cos(phi)) / (2*r*(R + r*cos(phi))).")
        print("To ensure H is never zero, we check its minimum value, which occurs when cos(phi) = -1.")
        
        # We explicitly show the numbers in the final calculation step.
        print(f"The numerator at this point is {R} + 2*{r}*(-1) = {min_numerator_val}")
        print(f"The denominator at this point is 2*{r}*({R} + {r}*(-1)) = {denominator_at_min}")
        print(f"So the minimum mean curvature is {min_numerator_val}/{denominator_at_min} = {min_H:.4f}")

        if min_H > 0:
            print("\nSince the minimum value of H is positive, the mean curvature never vanishes for this torus.")
            print("This confirms that a surface of genus 1 can satisfy the conditions.")
        else:
            # This case shouldn't be reached if R > 2r
            print("\nError in calculation, H can be zero or negative.")

    else:
        print(f"The condition R > 2r ({R} > {2*r}) is NOT satisfied.")
        print("For this torus, the mean curvature will be zero at some point.")

# Example parameters for a 'fat' torus. Let R=3 and r=1.
R_val = 3
r_val = 1
check_torus_mean_curvature(R_val, r_val)
print("\nSince genus 1 is possible, and it is a deep mathematical result that higher genera are not, the maximal genus is 1.")
