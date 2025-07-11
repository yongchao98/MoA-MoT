import math

def solve_representation_percentage():
    """
    Calculates the percentage of irreducible objects in the category C
    of indecomposable finite-dimensional u_q(sl_2) representations for q
    a primitive third root of unity.
    """
    
    # The order 'l' for q = exp(2*pi*i/3) is 3, since (q^2)^3 = 1 and 3 is the smallest such integer.
    l = 3

    # The category of representations decomposes into blocks based on the central character chi=(a,b,c).
    # We analyze the ratio of irreducible to indecomposable objects in each type of block.

    # Case 1: The principal block (a=0, b=0).
    # This block has a finite number of irreducible objects but infinitely many indecomposable ones.
    # The irreducible objects are L_0, L_1, L_2.
    num_irreducible_principal = 3
    num_indecomposable_principal = float('inf')
    # The percentage is effectively zero.
    perc_principal = 0.0

    # Case 2: Taft-like blocks (a=0, b!=0 or a!=0, b=0).
    # These blocks are equivalent to the representation category of the Taft algebra T_l.
    # T_l has l simple objects and l^2 indecomposable objects.
    num_irreducible_taft = l
    num_indecomposable_taft = l**2
    perc_taft = (num_irreducible_taft / num_indecomposable_taft) * 100

    # Case 3: Semi-simple blocks (a!=0, b!=0).
    # In these blocks, every indecomposable object is also irreducible.
    # Each such block contains a single indecomposable object, which is irreducible.
    num_irreducible_semisimple = 1
    num_indecomposable_semisimple = 1
    perc_semisimple = (num_irreducible_semisimple / num_indecomposable_semisimple) * 100
    
    print("The category C of indecomposable representations is partitioned into blocks based on a central character (a, b, c).")
    print("We analyze the percentage of irreducible objects within each type of block:\n")
    
    print(f"1. Principal Block (a=0, b=0):")
    print(f"   - Irreducible objects: {num_irreducible_principal}")
    print(f"   - Indecomposable objects: infinity")
    print(f"   - Percentage: {perc_principal:.1f}%\n")
    
    print(f"2. Taft-like Blocks (one of a, b is zero, but not both):")
    print(f"   - Irreducible objects: {num_irreducible_taft}")
    print(f"   - Indecomposable objects: {num_indecomposable_taft}")
    print(f"   - Percentage: ({num_irreducible_taft}/{num_indecomposable_taft}) * 100 = {perc_taft:.2f}%\n")

    print(f"3. Semi-simple Blocks (a!=0, b!=0):")
    print(f"   - Irreducible objects: {num_irreducible_semisimple}")
    print(f"   - Indecomposable objects: {num_indecomposable_semisimple}")
    print(f"   - Percentage: ({num_irreducible_semisimple}/{num_indecomposable_semisimple}) * 100 = {perc_semisimple:.1f}%\n")

    print("The space of all central characters (a, b, c) is C^3.")
    print("The characters corresponding to semi-simple blocks (where a!=0 and b!=0) form a Zariski-open dense subset.")
    print("Therefore, a 'generic' indecomposable representation belongs to a semi-simple block.")
    print("The overall percentage is thus determined by this generic case.\n")
    
    print("Final Equation:")
    print(f"Percentage = (Number of irreducibles in a generic block / Number of indecomposables in a generic block) * 100")
    print(f"Percentage = ({num_irreducible_semisimple} / {num_indecomposable_semisimple}) * 100 = {perc_semisimple:.0f}")

solve_representation_percentage()