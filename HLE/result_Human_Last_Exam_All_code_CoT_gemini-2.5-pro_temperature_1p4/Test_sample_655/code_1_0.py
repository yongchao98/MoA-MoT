import math

def solve():
    """
    Solves the energy ball packing problem by finding the optimal container and ball quantities.
    """
    max_sa = 1050.0
    
    # Step 3: Find the optimal integer dimensions (Nx, Ny, Nz) for packing 2-cm balls.
    # The goal is to maximize product Nx*Ny*Nz given Nx*Ny + Nx*Nz + Ny*Nz <= 32.8125
    
    max_balls = 0
    best_dims = (0, 0, 0)
    limit = 32.8125

    print("Step 1: Searching for the optimal number of 2-cm balls (Nx, Ny, Nz) to pack in a box.")
    print(f"The constraint is Nx*Ny + Nx*Nz + Ny*Nz <= {limit}\n")

    # We can assume Nx <= Ny <= Nz without loss of generality to avoid duplicate permutations.
    # Max possible value for Nx would be when Nx=Ny=Nz: 3*Nx^2 <= 32.8 -> Nx <= sqrt(10.9) -> Nx <= 3
    for nx in range(1, 4):
        # Similarly, 2*ny^2 <= 32.8 -> ny <= 4
        for ny in range(nx, 5):
            # from nx*ny + nz*(nx+ny) <= limit -> nz <= (limit - nx*ny) / (nx+ny)
            if (nx + ny) == 0: continue
            nz_limit = (limit - nx*ny) / (nx+ny)
            for nz in range(ny, int(nz_limit) + 2):
                if nz == 0: continue
                sa_term = nx * ny + nx * nz + ny * nz
                if sa_term <= limit:
                    num_balls = nx * ny * nz
                    if num_balls > max_balls:
                        max_balls = num_balls
                        best_dims = (nx, ny, nz)

    nx, ny, nz = best_dims
    num_2cm_balls = max_balls

    l = nx * 4
    w = ny * 4
    h = nz * 4
    
    surface_area = 2 * (l * w + l * h + w * h)
    
    print("Step 2: Optimal configuration for 2-cm balls found.")
    print(f"Ball counts along dimensions (Nx, Ny, Nz): ({nx}, {ny}, {nz})")
    print(f"Total 2-cm balls (b): {num_2cm_balls}")
    print(f"Resulting container: box {l}x{w}x{h}")
    print(f"Calculated surface area: {surface_area} cm^2 (must be <= 1050 cm^2)\n")

    # Step 4: Calculate the number of 1-cm balls that can fit in the interstitial spaces.
    # The number of such spaces is (Nx-1) * (Ny-1) * (Nz-1).
    num_1cm_balls = (nx - 1) * (ny - 1) * (nz - 1)
    
    # We need to check if these balls are valid (don't overlap with large balls or walls).
    # Small ball radius R1=1, large ball radius R2=2. Required center distance >= R1+R2=3.
    # The interstitial site is at the center of a cube of large balls.
    # Its distance to any of the 8 corner large-ball centers is sqrt(2^2+2^2+2^2) = sqrt(12) ~= 3.46.
    # Since 3.46 > 3, there is no overlap with large balls.
    # The walls are checked by ensuring the centers of small balls are > 1cm from them.
    # A center at (4,4,4) in a box starting at (0,0,0) is valid. This holds for our case.
    
    print("Step 3: Calculating how many 1-cm balls can fit in the gaps.")
    print("The largest gaps are at the center of the cubic cells formed by the 2-cm balls.")
    print(f"Number of such gaps: ({nx}-1) * ({ny}-1) * ({nz}-1) = {num_1cm_balls}")
    print(f"Total 1-cm balls (a): {num_1cm_balls}\n")
    
    # Step 5: Final result.
    total_energy = num_2cm_balls * 10 + num_1cm_balls * 1
    
    container_str = f"box {l}x{w}x{h}"
    num_1cm_str = f"{num_1cm_balls}"
    num_2cm_str = f"{num_2cm_balls}"
    
    print("Step 4: Final calculation and summary.")
    print("Container description (C): " + container_str)
    print("Number of 1-cm balls (a): " + num_1cm_str)
    print("Number of 2-cm balls (b): " + num_2cm_str)
    print(f"Total maximized energy: {num_2cm_str} * 10 MJ + {num_1cm_str} * 1 MJ = {total_energy} MJ")
    
    final_answer = f"[{container_str}]{num_1cm_str};{num_2cm_str}"
    print("\nFinal Answer Format: [C]a;b")
    print(final_answer)
    return final_answer

final_answer_string = solve()

# The final answer in the required format
print(f"\n<<<{final_answer_string}>>>")