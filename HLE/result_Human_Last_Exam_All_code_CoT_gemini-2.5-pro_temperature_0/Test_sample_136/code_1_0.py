import math

def solve_probability():
    """
    Calculates the limit of the conditional probability P[x_0 not visited | 0 not visited].

    The solution follows these steps:
    1. The limit is known to be equal to (P[LERW from 0 to infinity avoids x_0])^2.
    2. P[LERW avoids x_0] = 1 - P[LERW hits x_0].
    3. P[LERW hits x_0] = P_0[T_{x_0} < T_0^+] + P_{x_0}[T_0 < T_{x_0}^+].
       By symmetry, this is 2 * p, where p = P_0[T_{x_0} < T_0^+].
    4. p is calculated using committor probabilities h(y) for neighbors y of 0.
    5. Committor probabilities are calculated using known effective resistance values on the Z^2 grid.
    """

    # Step 1: Define problem parameters
    # 0 is at (0,0).
    # x_0 has two common neighbors with 0. This corresponds to x_0 = (1,1).
    # The neighbors of 0 are N(0) = {(1,0), (0,1), (-1,0), (0,-1)}.
    # The common neighbors of 0 and (1,1) are (1,0) and (0,1).
    
    print("Step 1: Define vertices and neighbors.")
    origin = (0, 0)
    x0 = (1, 1)
    neighbors_of_0 = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    print(f"Vertex 0 is at {origin}, vertex x_0 is at {x0}.")
    print(f"Neighbors of 0 are: {neighbors_of_0}\n")

    # Step 2: Use known effective resistance values R(u,v) on the infinite 2D grid.
    # These values are from the physics/mathematics literature (e.g., Wu, Glasser).
    # R(u,v) is the resistance between u and v. R(v1,v2) = R(v1-v2, 0).
    # We need R(1,0), R(1,1), R(2,1).
    R_1_0 = 0.5
    R_1_1 = 2 / math.pi
    R_2_1 = 2 / math.pi
    
    print("Step 2: Use known effective resistance values from Z^2 grid.")
    print(f"R((0,0), (1,0)) = {R_1_0}")
    print(f"R((0,0), (1,1)) = {R_1_1}")
    print(f"R((0,0), (2,1)) = {R_2_1}\n")

    # Step 3: Calculate committor probabilities h(y) = P_y[T_{x_0} < T_0] for y in N(0).
    # Formula: h(y) = 0.5 + (R(y,0) - R(y,x_0)) / (2 * R(0,x_0))
    
    print("Step 3: Calculate committor probabilities h(y) for each neighbor of 0.")
    R_0_x0 = R_1_1
    
    # For y = (1,0) (a common neighbor)
    # R((1,0), 0) = R(1,0)
    # R((1,0), x_0) = R((1,0), (1,1)) = R(0,1) = R(1,0)
    R_y_0 = R_1_0
    R_y_x0 = R_1_0
    h_1_0 = 0.5 + (R_y_0 - R_y_x0) / (2 * R_0_x0)
    print(f"For neighbor y=(1,0): h(y) = 0.5 + ({R_y_0} - {R_y_x0}) / (2 * {R_0_x0}) = {h_1_0}")

    # For y = (0,1) (a common neighbor)
    # R((0,1), 0) = R(0,1) = R(1,0)
    # R((0,1), x_0) = R((0,1), (1,1)) = R(1,0)
    R_y_0 = R_1_0
    R_y_x0 = R_1_0
    h_0_1 = 0.5 + (R_y_0 - R_y_x0) / (2 * R_0_x0)
    print(f"For neighbor y=(0,1): h(y) = 0.5 + ({R_y_0} - {R_y_x0}) / (2 * {R_0_x0}) = {h_0_1}")

    # For y = (-1,0)
    # R((-1,0), 0) = R(1,0)
    # R((-1,0), x_0) = R((-1,0), (1,1)) = R(2,1)
    R_y_0 = R_1_0
    R_y_x0 = R_2_1
    h_m1_0 = 0.5 + (R_y_0 - R_y_x0) / (2 * R_0_x0)
    print(f"For neighbor y=(-1,0): h(y) = 0.5 + ({R_y_0} - {R_y_x0}) / (2 * {R_0_x0}) = {h_m1_0}")

    # For y = (0,-1)
    # R((0,-1), 0) = R(0,1) = R(1,0)
    # R((0,-1), x_0) = R((0,-1), (1,1)) = R(1,2) = R(2,1)
    R_y_0 = R_1_0
    R_y_x0 = R_2_1
    h_0_m1 = 0.5 + (R_y_0 - R_y_x0) / (2 * R_0_x0)
    print(f"For neighbor y=(0,-1): h(y) = 0.5 + ({R_y_0} - {R_y_x0}) / (2 * {R_0_x0}) = {h_0_m1}\n")

    # Step 4: Calculate escape probability p = P_0[T_{x_0} < T_0^+]
    committor_probs = [h_1_0, h_0_1, h_m1_0, h_0_m1]
    p = sum(committor_probs) / 4
    print("Step 4: Calculate the escape probability p = (1/4) * sum(h(y)).")
    print(f"p = (1/4) * ({h_1_0} + {h_0_1} + {h_m1_0} + {h_0_m1}) = {p}\n")

    # Step 5: Calculate the LERW avoidance probability
    # P(LERW avoids x_0) = 1 - 2*p
    lerw_avoids_prob = 1 - 2 * p
    print("Step 5: Calculate the LERW avoidance probability = 1 - 2*p.")
    print(f"P(LERW avoids x_0) = 1 - 2 * {p} = {lerw_avoids_prob}\n")

    # Step 6: The final answer is the square of the LERW avoidance probability
    final_answer = lerw_avoids_prob ** 2
    print("Step 6: The final answer is the square of the LERW avoidance probability.")
    print(f"Final Limit = (P(LERW avoids x_0))^2")
    print(f"Final Limit = ({lerw_avoids_prob})^2 = {final_answer}")
    
    # Print the final equation with numbers
    print("\nThe final result is derived from the expression (1/2 - pi/8)^2:")
    val_1 = 1
    val_2 = 2
    val_pi = math.pi
    val_8 = 8
    result = (val_1/val_2 - val_pi/val_8)**2
    print(f"({val_1}/{val_2} - {val_pi}/{val_8})^2 = {result}")


solve_probability()