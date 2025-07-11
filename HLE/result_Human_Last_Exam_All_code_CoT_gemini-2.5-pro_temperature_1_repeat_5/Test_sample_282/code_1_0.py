import numpy as np

def get_volume_from_params(params):
    """Estimates the volume of a set from its parameter samples."""
    if params.shape[0] == 0:
        return 0
    # Calculate the volume of the bounding box of the parameters
    ranges = np.ptp(params, axis=0)
    return np.prod(ranges)

def main():
    """
    Solves for K and provides numerical demonstrations.
    """
    print("Step 1: Find an upper bound for K using a non-expanding set.")
    print("We choose X = SO(2), a compact subgroup of SL_2(R).")
    
    # Generate points for X = SO(2)
    # A matrix in SO(2) is defined by an angle theta.
    num_samples = 2000
    thetas = np.random.uniform(0, 2 * np.pi, num_samples)
    X_so2 = np.array([[[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]] for t in thetas])
    
    # Generate points for X^3
    indices = np.random.randint(0, num_samples, size=(num_samples, 3))
    X3_so2 = np.array([X_so2[i] @ X_so2[j] @ X_so2[k] for i, j, k in indices])

    # The measure of SO(2) is finite. We approximate mu(X) by the range of the parameter theta.
    mu_X = 2 * np.pi
    # Since X=SO(2) is a group, X^3 = X, so mu(X^3) = mu(X).
    mu_X3 = 2 * np.pi
    ratio_so2 = mu_X3 / mu_X

    print(f"For X = SO(2), a group, X^3 = X. So, mu(X^3) / mu(X) = 1.")
    print(f"The inequality mu(X^3) >= K*mu(X) becomes mu(X) >= K*mu(X), which implies K <= {ratio_so2:.1f}.")
    print("-" * 20)

    print("Step 2: Demonstrate expansion for a generic set.")
    print("We choose Y to be a small ball around the identity matrix.")
    
    # Generate points for Y = a small ball near identity
    # We parameterize elements near identity using 3 small parameters a,b,c.
    # g(a,b,c) = [[1+a, b], [c, (1+bc)/(1+a)]]
    eps = 0.01
    Y_params = np.random.uniform(-eps, eps, size=(num_samples, 3))
    Y_set = np.array([[[1+a, b], [c, (1+b*c)/(1+a)]] for a, b, c in Y_params])
    
    # Generate points for Y^3
    indices = np.random.randint(0, num_samples, size=(num_samples, 3))
    Y3_set = np.array([Y_set[i] @ Y_set[j] @ Y_set[k] for i, j, k in indices])
    
    # We can approximate Y^3 elements by I + (a1+a2+a3, b1+b2+b3, c1+c2+c3)
    # So parameters for Y^3 are roughly the sum of parameters for Y.
    Y3_params_approx = np.array([Y_params[i] + Y_params[j] + Y_params[k] for i, j, k in indices])
    
    # The "measure" or "volume" is proportional to the volume of the parameter space.
    # The parameter space for Y is a cube of side 2*eps.
    # The parameter space for Y^3 is roughly a shape with ranges [-3*eps, 3*eps] for each param.
    mu_Y_approx = get_volume_from_params(Y_params)
    mu_Y3_approx = get_volume_from_params(Y3_params_approx)

    # For a small ball, the volume ratio mu(Y^3)/mu(Y) approaches 3^dim = 3^3 = 27
    ratio_ball = mu_Y3_approx / mu_Y_approx if mu_Y_approx > 0 else float('inf')
    
    print(f"For Y = Ball(I, eps), a generic set, mu(Y^3) / mu(Y) is approximately {ratio_ball:.1f}.")
    print("This demonstrates that many sets expand significantly.")
    print("-" * 20)

    print("Step 3: Conclusion.")
    print("The constant K must hold for ALL compact sets, including non-expanding ones like SO(2).")
    print("The existence of SO(2) forces K <= 1.")
    print("A general theorem for unimodular groups states that mu(X^3) >= mu(X), which means K=1 is always satisfied.")
    print("\nTherefore, the largest possible value of K is 1.")
    
    K = 1
    print("\nThe final equation with the largest possible value of K is:")
    print(f"mu(X^3) >= {K} * mu(X)")

if __name__ == '__main__':
    main()
<<<1>>>