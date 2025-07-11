import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for K in the parliament design.
    """
    # Step 1: Define the constants from the problem description.
    r1 = 3.0  # Initial radius of the first row in meters.
    row_depth = 1.5  # Depth per row in meters.
    h_seated = 1.0  # Height of a seated person in meters.

    # Step 2: Formulate the visibility constraint.
    # The constraint is that the line of sight from any observer's head to the
    # speaker's feet must not be blocked by the head of anyone seated in between.
    # The worst-case for this constraint occurs with the blocker who is closest
    # to the speaker, i.e., the person in the second row.

    # Speaker's feet position (S_feet):
    # Located at radius r1, on the floor of the paraboloid.
    # (r, h) = (r1, r1^2 / K)
    
    # Blocker's head position (B_head):
    # The most restrictive blocker is in the second row (m=2).
    r_blocker = r1 + row_depth
    # Blocker's head is at radius r_blocker, 1m above the floor.
    # (r, h) = (r_blocker, r_blocker^2 / K + h_seated)

    print("To solve for K, we establish the critical visibility constraint.")
    print("An observer's line of sight to the speaker's feet must clear the head of any person in between.")
    print("This leads to the inequality: K < (r_blocker - r_speaker)^2")
    print("\n--- Calculation Steps ---")
    
    # Step 3: Set up and solve the inequality for K.
    # The derivation shows that for the view to be clear, K must satisfy:
    # K < (r_blocker - r1)^2
    
    print(f"The speaker is in the first row at a radius r_speaker = {r1} m.")
    print(f"The first potential blocker is in the second row at a radius r_blocker = {r1} + {row_depth} = {r_blocker} m.")
    
    # The difference in radius is the base of the triangle for the sightline calculation.
    delta_r = r_blocker - r1
    
    # The inequality is K < delta_r^2
    k_limit = delta_r ** 2
    
    print("\nThe derived inequality is K < (r_blocker - r_speaker)^2")
    print(f"Substituting the values: K < ({r_blocker} - {r1})^2")
    print(f"This simplifies to: K < ({delta_r})^2")
    print(f"So, the condition is: K < {k_limit}")
    
    # Step 4: Find the maximum integer value for K.
    # Since K must be an integer, the maximum value is the floor of k_limit.
    max_k = math.floor(k_limit)
    
    print(f"\nAs K must be an integer, the maximum value it can take is {max_k}.")
    print("\n--- Final Equation ---")
    print(f"K < ({r_blocker} - {r1})^2 = {k_limit}")
    print(f"Maximum Integer K = {max_k}")
    
    print(f"\n<<< {max_k} >>>")

solve_parliament_design()