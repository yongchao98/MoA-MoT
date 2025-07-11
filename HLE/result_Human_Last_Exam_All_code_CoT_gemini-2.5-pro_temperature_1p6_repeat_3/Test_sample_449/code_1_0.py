import math

def solve():
    """
    Calculates the approximate probability that a 2D random walk conditioned to avoid the origin,
    starting from (3000, 4000), will never hit the set of the origin's four neighbours.
    """
    # Step 1: Define the starting position and calculate its distance from the origin.
    x0, y0 = 3000, 4000
    r0 = math.sqrt(x0**2 + y0**2)

    # Step 2: Calculate the Newtonian capacity of the set A = {(1,0), (-1,0), (0,1), (0,-1)}.
    # The set is modeled as four points in 3D space: y1=(1,0,0), y2=(-1,0,0), y3=(0,1,0), y4=(0,-1,0).
    # The capacity C is given by C = 4q, where q is the charge on each point such that the potential is 1.
    # The potential at one point (e.g., y1) due to the others is:
    # V1 = q * (1/||y1-y2|| + 1/||y1-y3|| + 1/||y1-y4||)
    # ||y1-y2|| = 2
    # ||y1-y3|| = sqrt(2)
    # ||y1-y4|| = sqrt(2)
    # So, 1 = q * (1/2 + 1/sqrt(2) + 1/sqrt(2)) = q * (0.5 + sqrt(2)).
    # This leads to q = 1 / (0.5 + sqrt(2)).
    # The capacity C = 4q = 4 / (0.5 + sqrt(2)) = 8 / (1 + 2*sqrt(2)).
    sqrt_2 = math.sqrt(2)
    capacity = 8 / (1 + 2 * sqrt_2)

    # Step 3: Calculate the probability of hitting the set A.
    # P(hit) approx capacity / r0
    prob_hit = capacity / r0

    # Step 4: Calculate the probability of never hitting the set A.
    prob_never_hit = 1 - prob_hit

    # Print the equation and the final result.
    # The equation is: 1 - (Capacity / Distance)
    # where Capacity = 8 / (1 + 2 * sqrt(2))
    # and Distance = 5000
    print("The probability is calculated as: 1 - (C / r0)")
    print(f"where the distance r0 = {r0}")
    print(f"and the capacity C = 8 / (1 + 2 * sqrt(2)) ≈ {capacity:.4f}")
    
    equation = f"P(never hit) = 1 - ({capacity:.4f} / {r0}) = {prob_never_hit:.8f}"
    
    # We are asked to provide the result with two significant digits.
    # The value is 0.999582..., which rounds to 1.0 for two significant digits.
    final_answer = round(prob_never_hit, 2)
    if final_answer == 1.0:
        final_answer = "1.0"
        
    print(f"\nFinal Equation:")
    print(f"1 - ( {8} / (1 + 2 * {sqrt_2}) ) / {r0} = {prob_never_hit}")
    
    print(f"\nThe approximate probability rounded to two significant digits is: {final_answer}")
    
solve()
<<<1.0>>>