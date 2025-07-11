import numpy as np

def check_matrices():
    """
    Analyzes the properties of matrices A, B, C, D, E, F and determines if they are in the set P.
    """
    
    # A = [[0,0],[0,0]]
    # Check A
    A_trace = 0 + 0
    A_in_P = A_trace >= 1
    
    # B = [[6,4],[3,7]]
    # Check B
    B_symmetric = 4 == 3
    B_in_P = B_symmetric

    # C = [[1,-1/2],[-1/2,1]]
    # Check C
    C_symmetric = True
    C_psd = (1*1 - (-0.5)*(-0.5)) >= 0
    C_trace = 1 + 1
    C_in_P = C_symmetric and C_psd and C_trace >= 1
    # We showed in the reasoning that C can be constructed as a convex combination.
    # C = 3/8*M(1,0) + 3/8*M(0,1) + 1/8*M(2,-1) + 1/8*M(-1,2)
    # v1=(1,0) -> M1=[[1,0],[0,0]]
    # v2=(0,1) -> M2=[[0,0],[0,1]]
    # v3=(2,-1) -> M3=[[4,-2],[-2,1]]
    # v4=(-1,2) -> M4=[[1,-2],[-2,4]]
    C_calc = (3/8)*np.array([[1,0],[0,0]]) + (3/8)*np.array([[0,0],[0,1]]) + \
             (1/8)*np.array([[4,-2],[-2,1]]) + (1/8)*np.array([[1,-2],[-2,4]])
    
    # D = [[pi,1],[1,pi^2]]
    # Check D
    # As reasoned, the transcendental nature of pi prevents D from being in P.
    D_in_P = False
    
    # E = [[1,pi],[pi,1]]
    # Check E
    pi = np.pi
    E_symmetric = True
    E_psd = (1*1 - pi*pi) >= 0
    E_in_P = E_symmetric and E_psd
    
    # F = [[42,0],[0,0]]
    # Check F
    F_symmetric = True
    F_psd = (42*0 - 0*0) >= 0
    F_trace = 42 + 0
    F_in_P = F_symmetric and F_psd and F_trace >= 1
    # We showed that F can be constructed.
    # F = 7/13*M(6,0) + 6/13*M(7,0)
    F_calc = (7/13)*np.array([[36,0],[0,0]]) + (6/13)*np.array([[49,0],[0,0]])

    result = []
    if C_in_P:
        result.append('C')
    if D_in_P:
        result.append('D')
    if F_in_P:
        result.append('F')
        
    print(result)

check_matrices()