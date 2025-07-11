import numpy as np
import re

def solve_quantum_puzzle():
    """
    Solves the quantum transformation puzzle by checking for linearity.
    """
    # 1. Define the six basis states as complex vectors
    SQRT2 = np.sqrt(2)
    states = {
        '0': np.array([1, 0], dtype=complex),
        '1': np.array([0, 1], dtype=complex),
        '+': np.array([1/SQRT2, 1/SQRT2], dtype=complex),
        '-': np.array([1/SQRT2, -1/SQRT2], dtype=complex),
        'i': np.array([1/SQRT2, 1j/SQRT2], dtype=complex),
        '-i': np.array([1/SQRT2, -1j/SQRT2], dtype=complex),
    }

    # Coefficients to express superposition states in the |0>, |1> basis
    # e.g., |+> = (1/sqrt2)|0> + (1/sqrt2)|1>
    coeffs = {
        '+': np.array([1/SQRT2, 1/SQRT2], dtype=complex),
        '-': np.array([1/SQRT2, -1/SQRT2], dtype=complex),
        'i': np.array([1/SQRT2, 1j/SQRT2], dtype=complex),
        '-i': np.array([1/SQRT2, -1j/SQRT2], dtype=complex),
    }

    # The list of transformations for each answer choice
    transforms_data = """
    A. |0> -> |1>   |1> -> |0>   |+> -> |->   |-> -> |+>   |i> -> |i>   |-i> -> |-i>
    B. |0> -> |0>   |1> -> |1>   |+> -> |->   |-> -> |+>   |i> -> |-i>  |-i> -> |i>
    C. |0> -> |+>   |1> -> |->   |+> -> |0>   |-> -> |1>   |i> -> |-i>  |-i> -> |i>
    D. |0> -> |0>   |1> -> |1>   |+> -> |i>   |-> -> |-i>  |i> -> |->   |-i> -> |+>
    E. |0> -> |0>   |1> -> |1>   |+> -> |+>   |-> -> |->   |i> -> |i>   |-i> -> |-i>
    F. |0> -> |+>   |1> -> |->   |+> -> |1>   |-> -> |0>   |i> -> |i>   |-i> -> |-i>
    G. |0> -> |->   |1> -> |+>   |+> -> |0>   |-> -> |1>   |i> -> |i>   |-i> -> |-i>
    H. |0> -> |i>   |1> -> |-i>  |+> -> |+>   |-> -> |->   |i> -> |1>   |-i> -> |0>
    I. |0> -> |-i>  |1> -> |i>   |+> -> |+>   |-> -> |->   |i> -> |0>   |-i> -> |1>
    J. |0> -> |+>   |1> -> |->   |+> -> |-i>  |-> -> |i>   |i> -> |1>   |-i> -> |0>
    K. |0> -> |-i>  |1> -> |i>   |+> -> |0>   |-> -> |1>   |i> -> |->   |-i> -> |+>
    L. |0> -> |i>   |1> -> |-i>  |+> -> |0>   |-> -> |1>   |i> -> |+>   |-i> -> |->
    M. |0> -> |+>   |1> -> |->   |+> -> |i>   |-> -> |-i>  |i> -> |0>   |-i> -> |1>
    N. |0> -> |->   |1> -> |+>   |+> -> |i>   |-> -> |-i>  |i> -> |1>   |-i> -> |0>
    O. |0> -> |-i>  |1> -> |i>   |+> -> |1>   |-> -> |0>   |i> -> |+>   |-i> -> |->
    P. |0> -> |->   |1> -> |+>   |+> -> |-i>  |-> -> |i>   |i> -> |0>   |-i> -> |1>
    Q. |0> -> |->   |1> -> |+>   |+> -> |-i>  |-> -> |i>   |i> -> |1>   |-i> -> |0>
    R. |0> -> |i>   |1> -> |-i>  |+> -> |1>   |-> -> |0>   |i> -> |->   |-i> -> |+>
    S. |0> -> |i>   |1> -> |-i>  |+> -> |->   |-> -> |+>   |i> -> |0>   |-i> -> |1>
    T. |0> -> |-i>  |1> -> |i>   |+> -> |->   |-> -> |+>   |i> -> |1>   |-i> -> |0>
    U. |0> -> |->   |1> -> |+>   |+> -> |1>   |-> -> |0>   |i> -> |-i>  |-i> -> |i>
    V. |0> -> |1>   |1> -> |0>   |+> -> |i>   |-> -> |-i>  |i> -> |+>   |-i> -> |->
    W. |0> -> |1>   |1> -> |0>   |+> -> |-i>  |-> -> |i>   |i> -> |->   |-i> -> |+>
    """
    
    def are_proportional(v1, v2, tol=1e-9):
        # Two 2D vectors v1, v2 are proportional if the determinant of the
        # matrix [v1, v2] is zero. This is more numerically stable than division.
        # det = v1[0]*v2[1] - v1[1]*v2[0]
        return abs(v1[0] * v2[1] - v1[1] * v2[0]) < tol

    all_transforms = {}
    for line in transforms_data.strip().split('\n'):
        line = line.strip()
        choice = line[0]
        # Find all transformations like '|in> -> |out>'
        pairs = re.findall(r'\|(.*?)\⟩\s*->\s*\|(.*?)\⟩', line)
        all_transforms[choice] = {in_state: out_state for in_state, out_state in pairs}

    # 2. Loop through each answer choice and check for linearity
    for choice, trans_map in all_transforms.items():
        is_possible = True
        
        # Get the output vectors from the transformation map
        out0_vec = states[trans_map['0']]
        out1_vec = states[trans_map['1']]
        out_plus_vec = states[trans_map['+']]

        # 3. Determine the required relative phase p_ratio = p1/p0
        # from the |+> transformation.
        # We need alpha*out0 + beta*p_ratio*out1 to be proportional to out_plus
        # alpha*out0[0] + beta*p_ratio*out1[0] = k * out_plus[0]
        # alpha*out0[1] + beta*p_ratio*out1[1] = k * out_plus[1]
        # Solving for p_ratio yields:
        alpha, beta = coeffs['+']
        num = alpha * (out0_vec[1] * out_plus_vec[0] - out0_vec[0] * out_plus_vec[1])
        den = beta * (out1_vec[0] * out_plus_vec[1] - out1_vec[1] * out_plus_vec[0])

        if abs(den) < 1e-9:
             # This case means out1 is proportional to out_plus.
             # For linearity, out0 must also be proportional.
             # But out0 and out1 are orthogonal for all given options, so this path is not taken.
             continue
        
        p_ratio = num / den
        
        # The relative phase must have a magnitude of 1
        if not np.isclose(abs(p_ratio), 1.0):
            is_possible = False
            failure_reason = "Derived relative phase is not a valid phase factor."
        else:
            # 4. Check consistency for the other superposition states
            for super_state_name in ['-', 'i', '-i']:
                alpha, beta = coeffs[super_state_name]
                
                # Calculate the output vector based on linearity
                calculated_out_vec = alpha * out0_vec + beta * p_ratio * out1_vec
                # Get the required output vector from the problem
                required_out_vec = states[trans_map[super_state_name]]
                
                if not are_proportional(calculated_out_vec, required_out_vec):
                    is_possible = False
                    print(f"The impossible transformation is: {choice}")
                    print(f"Reason: Linearity is violated for the input state |{super_state_name}⟩.")
                    print(f"Based on the transformations of |0> and |1>, the output for |{super_state_name}⟩ should be proportional to {np.round(calculated_out_vec, 3)}.")
                    print(f"However, the required output is proportional to {np.round(required_out_vec, 3)}.")
                    return # Exit after finding the first impossible transformation
    
solve_quantum_puzzle()