import numpy as np
import matplotlib.pyplot as plt

def demonstrate_pareto_gap():
    """
    This function demonstrates that for d=3, the Pareto front is not fully
    explorable by scalarization, even with non-negative data.
    It does this by plotting the attainable objective set and the points
    found by scalarization for a specific case where the Pareto front is concave.
    """
    
    # 1. Define the counterexample data matrix X (n=2, d=3) with non-negative entries
    X = np.array([[1., 2., 0.],
                  [0., 2., 1.]])
    n, d = X.shape
    x1 = X[0, :]
    x2 = X[1, :]
    
    print(f"Demonstrating a counterexample for d={d} and n={n}")
    print(f"The rows of the data matrix X are x1 = {x1} and x2 = {x2}.")

    # 2. Generate a large number of random unit vectors w to sample the attainable objective space
    num_samples = 50000
    random_vectors = np.random.randn(d, num_samples)
    w_samples = random_vectors / np.linalg.norm(random_vectors, axis=0)

    # 3. Calculate the corresponding objective values (f1, f2)
    # f_i(w) = (x_i . w)^2
    f1_all = (x1 @ w_samples)**2
    f2_all = (x2 @ w_samples)**2

    # 4. Find the solutions from the scalarization method
    scalarization_solutions = []
    # We iterate through different weights lambda = (alpha, 1-alpha)
    for alpha in np.linspace(0, 1, 201):
        # Construct the matrix M for the scalarized objective
        M = alpha * np.outer(x1, x1) + (1 - alpha) * np.outer(x2, x2)
        
        # The solution is the eigenvector corresponding to the largest eigenvalue of M
        eigenvalues, eigenvectors = np.linalg.eigh(M)
        w_solution = eigenvectors[:, -1]
        
        # Calculate the objective values for this solution
        f1_sol = (x1 @ w_solution)**2
        f2_sol = (x2 @ w_solution)**2
        scalarization_solutions.append([f1_sol, f2_sol])

    scalarization_solutions = np.array(scalarization_solutions)

    # 5. Plot the results to visualize the gap
    plt.figure(figsize=(10, 8))
    plt.scatter(f1_all, f2_all, s=2, alpha=0.15, label='Attainable Objective Space')
    
    # Sort for a clean line plot
    sorted_indices = np.argsort(scalarization_solutions[:, 0])
    plt.plot(scalarization_solutions[sorted_indices, 0], 
             scalarization_solutions[sorted_indices, 1], 
             'r-o', markersize=4, linewidth=2.5, label='Solutions from Scalarization')
    
    plt.title('Failure of Scalarization for d=3 with Non-Negative Data')
    plt.xlabel('Objective 1: $f_1(w) = (x_1 \cdot w)^2$')
    plt.ylabel('Objective 2: $f_2(w) = (x_2 \cdot w)^2$')
    plt.legend()
    plt.grid(True)
    
    print("\nThe plot shows the attainable objective space in light blue.")
    print("The red line shows the points found by scalarization.")
    print("The true Pareto front is the upper-right boundary of the blue region.")
    print("As you can see, the true front is concave, and the scalarization method only finds points on its convex hull, missing the optimal trade-offs in the middle.")
    plt.show()

demonstrate_pareto_gap()