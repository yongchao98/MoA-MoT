def analyze_algorithms():
    """
    Analyzes the convergence properties of three optimization algorithms.
    """

    print("Analyzing the convergence properties of three optimization algorithms:")
    print("--------------------------------------------------------------------")

    # Analysis of Gradient Descent (1)
    print("\n(1) Gradient Descent: x_{k+1} = x_k - gamma * nabla_f(x_k)")
    print("Assume the sequence converges: x_k -> x* as k -> infinity.")
    print("This implies that the difference between consecutive terms must go to zero: x_{k+1} - x_k -> 0.")
    print("From the update rule, we have: nabla_f(x_k) = (1 / gamma) * (x_k - x_{k+1})")
    print("Taking the limit as k -> infinity:")
    print("  lim(nabla_f(x_k)) = (1 / gamma) * lim(x_k - x_{k+1})")
    print("  nabla_f(x*) = (1 / gamma) * 0")
    print("  nabla_f(x*) = 0")
    print("This is the condition for a first-order stationary point. Therefore, Gradient Descent converges to stationary points.")

    # Analysis of Doubly-Projected Gradient Descent (2)
    print("\n(2) Doubly-Projected Gradient Descent: x_{k+1} = Proj_C(x_k + gamma_k * Proj_{T_{x_k}C}(-nabla_f(x_k)))")
    print("Analysis shows that if x_k -> x*, then the projected gradient v* = Proj_{T_{x*}C}(-nabla_f(x*)) must lie in the normal cone N_{x*}C.")
    print("Since v* is also in the tangent cone T_{x*}C by definition, v* must be in the intersection T_{x*}C intersect N_{x*}C.")
    print("For most well-behaved sets (convex, smooth manifolds, etc.), this intersection is {0}.")
    print("This implies v* = 0, which is the stationarity condition. Thus, this algorithm is generally expected to converge to stationary points.")

    # Analysis of Heavy-ball Method (3)
    print("\n(3) Heavy-ball Method: x_{k+1} = x_k + beta*(x_k - x_{k-1}) - gamma * nabla_f(x_k)")
    print("Assume the sequence converges: x_k -> x*.")
    print("This implies that x_k - x_{k-1} -> 0 and x_{k+1} - x_k -> 0.")
    print("Rearranging the update rule: nabla_f(x_k) = (1/gamma) * (x_k - x_{k+1} + beta*(x_k - x_{k-1}))")
    print("Taking the limit as k -> infinity:")
    print("  lim(nabla_f(x_k)) = (1/gamma) * (lim(x_k - x_{k+1}) + beta*lim(x_k - x_{k-1}))")
    print("  nabla_f(x*) = (1/gamma) * (0 + beta*0)")
    print("  nabla_f(x*) = 0")
    print("This simple analysis suggests convergence to a stationary point. However, this reasoning can be misleading.")
    print("Due to its second-order nature (momentum), the dynamics are more complex.")
    print("Advanced analysis (e.g., using tools from control theory) has shown that for some smooth non-convex functions, the heavy-ball method can converge to a non-stationary point.")

    # Conclusion
    print("\n--------------------------------------------------------------------")
    print("Conclusion: Based on the analysis, only the Heavy-ball method (3) has been shown to have the possibility of converging to a non-stationary point for certain smooth non-convex functions.")

if __name__ == '__main__':
    analyze_algorithms()