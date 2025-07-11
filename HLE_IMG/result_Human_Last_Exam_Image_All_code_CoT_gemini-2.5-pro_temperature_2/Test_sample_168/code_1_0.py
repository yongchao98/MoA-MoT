import collections
from itertools import permutations

def get_mirror_type(i, j, step):
    """Determines the mirror type for a step from cell (i,j)."""
    if step == 'L':
        # Crossing the left wall of cell (i,j)
        return 'G1' if i % 2 == 0 else 'G3'
    elif step == 'R':
        # Crossing the right wall of cell (i,j)
        return 'G3' if i % 2 == 0 else 'G1'
    elif step == 'D':
        # Crossing the bottom wall of cell (i,j)
        return 'G2' if j % 2 == 0 else 'G4'
    elif step == 'U':
        # Crossing the top wall of cell (i,j)
        return 'G4' if j % 2 == 0 else 'G2'
    return None

def solve():
    """
    Calculates the number of ways a light ray can reflect according to the problem's conditions.
    """
    # The set of steps corresponding to the reflections.
    # 2 on G1(L), 1 on G3(R), 1 on G2(D), 1 on G4(U)
    steps_to_permute = ['L', 'L', 'R', 'U', 'D']
    
    # The target multiset of mirrors.
    target_mirrors = collections.Counter({'G1': 2, 'G3': 1, 'G2': 1, 'G4': 1})
    
    valid_paths_count = 0
    
    # Iterate through all unique permutations of the steps.
    unique_permutations = set(permutations(steps_to_permute))
    
    for path in unique_permutations:
        i, j = 0, 0
        path_mirrors = []
        
        # For the current path, trace the steps and determine the mirror sequence.
        for step in path:
            mirror = get_mirror_type(i, j, step)
            path_mirrors.append(mirror)
            
            # Update cell coordinates based on the step.
            if step == 'L':
                i -= 1
            elif step == 'R':
                i += 1
            elif step == 'D':
                j -= 1
            elif step == 'U':
                j += 1
                
        # Check if the generated mirror sequence matches the target.
        if collections.Counter(path_mirrors) == target_mirrors:
            valid_paths_count += 1
            
    print(f"The number of ways to draw the light ray path is: {valid_paths_count}")

solve()
<<<16>>>