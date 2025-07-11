import sys

def solve_topological_classification():
    """
    Solves for the topological invariant group of a defect in a free fermion system
    based on the tenfold way classification.
    """
    print("### Solving for the Topological Invariant Group ###")
    print("\nStep 1: Identify the symmetry class.")
    time_reversal_sq = -1
    particle_hole_sq = -1
    # Based on the Altland-Zirnbauer classification:
    # T^2=-1, P^2=-1 corresponds to Class DIII.
    symmetry_class = "DIII"
    print(f"A system with time-reversal symmetry T^2 = {time_reversal_sq} "
          f"and particle-hole symmetry P^2 = {particle_hole_sq} belongs to Class {symmetry_class}.")

    print("\nStep 2: Determine the defect codimension (D).")
    system_dimension = 2
    # The user specified 'point defect (codimension D=1)'. This is ambiguous.
    # A point defect in 2D has dimension 0, so its codimension D = 2 - 0 = 2.
    # A line defect in 2D has dimension 1, so its codimension D = 2 - 1 = 1.
    # We will assume 'point defect' was the primary intent, which is physically significant.
    defect_type = "point defect"
    codimension_D = 2
    print(f"The system is 2D and contains a {defect_type}.")
    print(f"The codimension of a point defect in a 2D space is D = (system dimension) - (defect dimension) = 2 - 0 = {codimension_D}.")
    print("We will proceed with D=2 based on the term 'point defect'.")

    print("\nStep 3: Determine the effective dimension for classification (d_eff).")
    # The classification of a defect of codimension D is given by the bulk classification
    # in d_eff = D - 1 dimensions for the same symmetry class.
    effective_dimension = codimension_D - 1
    print("The rule for defect classification maps it to a bulk classification in d_eff = D - 1 dimensions.")
    print(f"Calculation: d_eff = {codimension_D} - 1 = {effective_dimension}")

    print(f"\nStep 4: Find the topological invariant for Class {symmetry_class} in d = {effective_dimension}.")
    # Periodic table for Real symmetry classes (condensed matter convention)
    # The key is the class name, the value is a list of invariants for d = 0, 1, 2, ...
    periodic_table = {
        "AI":    ["Z", "0", "0", "0", "Z2", "Z2", "0", "Z"],
        "BDI":   ["Z2", "Z", "0", "0", "0", "Z2", "Z2", "0"],
        "D":     ["Z2", "Z2", "Z", "0", "0", "0", "Z2", "Z2"],
        "DIII":  ["0", "Z2", "Z2", "Z", "0", "0", "0", "Z2"],
        "AII":   ["Z2", "0", "Z2", "Z2", "Z", "0", "0", "0"],
        "CII":   ["0", "Z2", "0", "Z2", "Z2", "Z", "0", "0"],
        "C":     ["0", "0", "Z2", "0", "Z2", "Z2", "Z", "0"],
        "CI":    ["0", "0", "0", "Z2", "0", "Z2", "Z2", "Z"]
    }
    # For pretty printing on one line, use unicode subscripts.
    sub = str.maketrans("012", "₀₁₂")
    if sys.stdout.encoding == 'utf-8':
        class_name_map = {"Z":"\u2124", "Z2":"\u2124\u2082"}
        final_group_str = class_name_map.get(
            periodic_table[symmetry_class][effective_dimension],
            periodic_table[symmetry_class][effective_dimension]
        )
    else: # Fallback for non-utf8 terminals
        final_group_str = periodic_table[symmetry_class][effective_dimension].replace("Z", "Z_")

    print(f"Looking up the entry for Class {symmetry_class} at d={effective_dimension} in the periodic table.")
    invariant_group = periodic_table[symmetry_class][effective_dimension]
    print(f"\nFinal Answer: The group of the topological invariant is {final_group_str}.")
    
    return invariant_group

# Run the solver and extract the final answer.
final_answer = solve_topological_classification()

# The final answer format for the platform
# It should be Z_2 or Z2. Using Z_2 for clarity
print(f"\n<<<{final_answer.replace('Z', 'Z_')}>>>")