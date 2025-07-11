def print_minimum_curvature_cost():
    """
    This function prints the equation for the minimum achievable curvature cost
    for the NGD update rule under the specified conditions.
    """
    
    # The minimum cost is achieved by using the Woodbury identity in combination with
    # an iterative solver like Conjugate Gradient (CG). The cost is measured in
    # floating-point operations (FLOPs).
    
    # The calculation involves matrix-vector products with the Jacobian J and its
    # transpose, avoiding the explicit formation of large matrices.
    # The total cost can be approximated by the following formula.
    
    # Let's break down the FLOPs:
    # 1. Compute J @ g:               ~2 * n * d^2
    # 2. k iterations of CG to solve (JJ^T + alpha*I)y = Jg:
    #    Each iteration needs J@(J^T@w), costing ~4 * n * d^2
    #    Total for k iterations:        ~4 * k * n * d^2
    # 3. Compute J^T @ y:               ~2 * n * d^2
    # 4. Final scaling and addition:    ~2 * d^2
    #
    # Total cost = (2*n*d^2) + (4*k*n*d^2) + (2*n*d^2) + (2*d^2)
    #            = (4*n*d^2 + 4*k*n*d^2) + 2*d^2
    #            = (4*k + 4)*n*d^2 + 2*d^2
    
    print("The minimum achievable curvature cost, derived from an efficient iterative algorithm, is given by the following equation in terms of floating-point operations (FLOPs):")
    
    # We output each number (coefficient) in the final equation.
    c1 = 4
    c2 = 4
    c3 = 2
    
    print(f"\nCost = ({c1}*k + {c2}) * n * d^2 + {c3}*d^2\n")
    
    print("Where:")
    print("  d: the size of the neural network layer.")
    print("  n: the number of training samples (with the condition n < d).")
    print("  k: the number of iterations for the Conjugate Gradient (CG) solver, which is typically a small number (k << n).")

if __name__ == '__main__':
    print_minimum_curvature_cost()
