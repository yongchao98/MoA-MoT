import numpy as np

def solve_inversion_barrier():
    """
    This function solves for the unknown inversion barrier by finding a
    quadratic relationship based on the provided data.
    """
    
    # N represents the number of five-membered rings in the molecules.
    # B is the corresponding inversion barrier in kcal/mol.
    
    # Data points from the problem description:
    # N=1 -> B=10
    # N=2 -> B=49
    # We add a physically sound baseline: N=0 (planar molecule) -> B=0.
    N_data = np.array([0, 1, 2])
    B_data = np.array([0, 10, 49])
    
    # We want to find a function B(N) = a*N^2 + b*N + c that fits the data.
    # np.polyfit finds the coefficients [a, b, c] for a polynomial of degree 2.
    coeffs = np.polyfit(N_data, B_data, 2)
    a, b, c = coeffs
    
    # The number of five-membered rings for the target molecule.
    N_target = 3
    
    # Calculate the predicted barrier using the derived formula.
    B_predicted = a * N_target**2 + b * N_target + c

    print("Step 1: Identified the number of five-membered rings (N) as the key structural variable.")
    print(f"    - Molecule 1 (Corannulene): N=1, Barrier={B_data[1]} kcal/mol")
    print(f"    - Molecule 2 (Diacenaphtho...): N=2, Barrier={B_data[2]} kcal/mol")
    print(f"    - Molecule 3 (Triacenaphtho...): N=3, Barrier=?")
    print("\nStep 2: Established a physical baseline: for N=0, the barrier B is 0 kcal/mol.")
    
    print("\nStep 3: Found a quadratic relationship B(N) = a*N^2 + b*N + c that fits the points (0,0), (1,10), and (2,49).")
    print(f"    - The calculated coefficients are: a = {a:.1f}, b = {b:.1f}, c = {c:.1f}")
    
    print("\nStep 4: Using the derived formula B(N) = {:.1f}*N^2 + {:.1f}*N, we can predict the barrier for N=3.".format(a,b))
    
    print("\nFinal Calculation:")
    print(f"Inversion Barrier = {a:.1f} * ({N_target})^2 + ({b:.1f}) * {N_target}")
    print(f"Inversion Barrier = {a:.1f} * {N_target**2} - {abs(b):.1f} * {N_target}")
    print(f"Inversion Barrier = {a * N_target**2:.1f} - {abs(b * N_target):.1f}")
    print(f"Inversion Barrier = {int(round(B_predicted))} kcal/mol")
    
    # The final answer in the requested format
    print(f"\n<<< {int(round(B_predicted))} >>>")

solve_inversion_barrier()