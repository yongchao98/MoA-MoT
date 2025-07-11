import math

def solve_pioneer_probe_packing():
    """
    This script calculates the optimal packing of energy balls for the Pioneer probe
    based on the analysis that a cubic box is the most efficient container shape.
    """

    # 1. Define the optimal container based on analysis
    shape_desc = "box"
    # Dimensions are chosen to maximize volume while staying under the surface area limit
    l, w, h = 13.0, 13.0, 13.0
    container_str = f"{shape_desc} {l}x{w}x{h}"

    surface_area = 2 * (l*w + w*h + h*l)

    # 2. Calculate the number of 2-cm radius balls (n2)
    # A 2-cm radius ball (4-cm diameter) is packed on a simple cubic lattice with 4cm spacing.
    # For a box of size L, center must be within [-L/2 + r, L/2 - r].
    # Box extent: [-6.5, 6.5]. Ball radius: 2.0. Center must be in [-4.5, 4.5].
    # Possible center coordinates (multiples of 0.5) along one axis with spacing 4.0:
    # -4.5, -0.5, 3.5  (3 positions)
    positions_per_axis_n2 = 3
    num_2cm_balls = positions_per_axis_n2 ** 3

    # 3. Calculate the number of 1-cm radius balls (n1)
    # The 3x3x3 lattice of 2-cm balls creates a 2x2x2 grid of empty cubes (voids).
    # We can place one 1-cm ball in the center of each of these voids.
    # The logic is verified because the distance from a void's center to the nearest
    # 2-cm ball center is sqrt(2^2 + 2^2 + 2^2) = ~3.46 cm, which is greater than the
    # required minimum distance of 2cm + 1cm = 3cm.
    positions_per_axis_n1 = positions_per_axis_n2 - 1
    num_1cm_balls = positions_per_axis_n1 ** 3

    # 4. Calculate the total energy
    energy_per_1cm_ball = 1
    energy_per_2cm_ball = 10
    total_energy = (num_1cm_balls * energy_per_1cm_ball) + \
                     (num_2cm_balls * energy_per_2cm_ball)

    # 5. Print the results and the final answer in the required format
    print(f"Optimal container found: A {container_str}")
    print(f"Surface Area: {surface_area:.2f} cm2 (Constraint: <= 1050 cm2)")
    print("-" * 30)
    print("Packing result:")
    print(f"Number of 2-cm radius balls (b): {num_2cm_balls}")
    print(f"Number of 1-cm radius balls (a): {num_1cm_balls}")
    print("-" * 30)
    print("Total potential energy equation:")
    # Here we output each number in the final equation, as requested.
    print(f"{num_2cm_balls} * {energy_per_2cm_ball} MJ + {num_1cm_balls} * {energy_per_1cm_ball} MJ = {total_energy} MJ")

    final_answer = f"[{container_str}]{num_1cm_balls};{num_2cm_balls}"
    print(f"\nFinal Answer: {final_answer}")
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    solve_pioneer_probe_packing()
