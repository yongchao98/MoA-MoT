import math

def solve_cube_cutting():
    """
    Calculates the minimum number of cuts to slice a 4x4x4 cube into 1x1x1 cubes
    with a knife of depth 2cm.
    """
    cube_size = 4
    knife_depth = 2
    
    # The positions of the required cut planes for each dimension
    cut_planes = list(range(1, cube_size)) # [1, 2, 3]
    
    center_cut_pos = cube_size / 2
    
    total_cuts = 0
    cut_costs = []
    
    print("Step-by-step calculation of the minimum number of cuts:\n")
    print(f"To cut a {cube_size}x{cube_size}x{cube_size} cube, we need {cube_size - 1} cut planes along each of the 3 axes.")
    print(f"The knife can only cut {knife_depth}cm deep. The cost of a cut is ceil(Thickness / {knife_depth}).\n")
    
    # Iterate through each of the 3 dimensions (X, Y, Z)
    for dim in ['X', 'Y', 'Z']:
        # Iterate through each required cut plane for that dimension
        for plane_pos in cut_planes:
            
            # Determine the thickness of the material for this cut
            if plane_pos == center_cut_pos:
                # Center cuts must go through the entire 4cm dimension
                thickness = cube_size
            else:
                # Outer cuts are made on pieces that are only 2cm thick
                thickness = cube_size / 2
            
            # Calculate the cost (number of knife actions) for this single plane
            cost = math.ceil(thickness / knife_depth)
            
            print(f"Cost for cut plane at {dim}={plane_pos}:")
            print(f"  - Thickness of material = {int(thickness)}cm")
            print(f"  - Number of actions = ceil({int(thickness)} / {knife_depth}) = {cost}")
            
            total_cuts += cost
            cut_costs.append(str(cost))

    print("\n-----------------------------------------")
    print("Final Calculation:")
    equation = " + ".join(cut_costs)
    print(f"Total Cuts = {equation} = {total_cuts}")
    print("-----------------------------------------")

solve_cube_cutting()
<<<12>>>