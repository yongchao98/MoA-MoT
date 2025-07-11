import numpy as np

def final_analysis():
    print("This script verifies that matrices C and F can be constructed as convex combinations of the generator matrices.")
    print("="*40)
    
    # Matrix C analysis
    print("For matrix C = [[1.0, -0.5], [-0.5, 1.0]]:")
    
    # We propose C = lambda1 * M_v1 + lambda2 * M_v2
    # With v1=(1,-1), v2=(1,1). We solve the system:
    # lambda1 + lambda2 = 1        (from diagonal entries)
    # -lambda1 + lambda2 = -0.5    (from off-diagonal entries)
    # Adding the two equations: 2*lambda2 = 0.5  => lambda2 = 0.25
    lambda2_c = 0.25
    lambda1_c = 1 - lambda2_c
    
    v1_c = np.array([[1], [-1]])
    v2_c = np.array([[1], [1]])
    M_v1_c = v1_c @ v1_c.T
    M_v2_c = v2_c @ v2_c.T
    C_reconstructed = lambda1_c * M_v1_c + lambda2_c * M_v2_c

    print(f"We construct C with lambda1 = {lambda1_c} and lambda2 = {lambda2_c}.")
    print(f"C = {lambda1_c} * (v1*v1^T) + {lambda2_c} * (v2*v2^T) where v1=[1,-1], v2=[1,1]")
    print(f"Equation: {lambda1_c} * [[{M_v1_c[0,0]}, {M_v1_c[0,1]}], [{M_v1_c[1,0]}, {M_v1_c[1,1]}]] + {lambda2_c} * [[{M_v2_c[0,0]}, {M_v2_c[0,1]}], [{M_v2_c[1,0]}, {M_v2_c[1,1]}]] = {C_reconstructed.tolist()}")
    print("The construction is valid, so C is in P.")

    print("="*40)

    # Matrix F analysis
    print("For matrix F = [[42, 0], [0, 0]]:")
    # For F's second diagonal element to be 0, all generating vectors must have b=0.
    # So we need to form 42 as a convex combination of squares a_i^2.
    # Let's use a_1^2 = 36 and a_2^2 = 49.
    # 42 = lambda1 * 36 + lambda2 * 49
    # lambda1 + lambda2 = 1 => lambda2 = 1 - lambda1
    # 42 = 36*lambda1 + 49*(1-lambda1) = 49 - 13*lambda1 => 13*lambda1 = 7
    lambda1_f = 7/13
    lambda2_f = 1 - lambda1_f

    v1_f = np.array([[6], [0]])
    v2_f = np.array([[7], [0]])
    M_v1_f = v1_f @ v1_f.T
    M_v2_f = v2_f @ v2_f.T
    F_reconstructed = lambda1_f * M_v1_f + lambda2_f * M_v2_f

    print(f"We construct F with lambda1 = {lambda1_f:.4f} and lambda2 = {lambda2_f:.4f}.")
    print(f"F = {lambda1_f:.4f} * (v1*v1^T) + {lambda2_f:.4f} * (v2*v2^T) where v1=[6,0], v2=[7,0]")
    print(f"{lambda1_f:.4f} * [[{M_v1_f[0,0]}, {M_v1_f[0,1]}], [{M_v1_f[1,0]}, {M_v1_f[1,1]}]] + {lambda2_f:.4f} * [[{M_v2_f[0,0]}, {M_v2_f[0,1]}], [{M_v2_f[1,0]}, {M_v2_f[1,1]}]] = {F_reconstructed.tolist()}")
    print("The construction is valid, so F is in P.")

    print("="*40)
    final_list = ['C', 'F']
    print(f"The list of matrices in P is: {final_list}")

final_analysis()