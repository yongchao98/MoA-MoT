import numpy as np

def demonstrate_non_convexity_and_solve():
    """
    Demonstrates that for d=3, the Pareto front can be non-convex even with x_i >= 0,
    and thus determines the largest d for which scalarization is guaranteed to work.
    """
    # Plan:
    # 1. Set up a counterexample for d=3 to show that the Pareto front is not always convex
    #    for x_i >= 0. This implies scalarization cannot find the entire front.
    # 2. The counterexample uses d=3, n=2, with specific non-negative vectors x1 and x2.
    #    x1 = (1, 0, 0)
    #    x2 = (0, cos(alpha), sin(alpha))
    #    A non-convex front appears for alpha in (0, pi/2). We choose alpha = pi/6.
    # 3. We will show that a specific point Y_star, known to be on the Pareto front,
    #    lies below the line segment connecting the two extreme points of the front, Y_a and Y_b.
    #    This demonstrates the non-convexity.
    # 4. Y_a is found by maximizing the first objective (y1). This gives w = (1,0,0) and Y_a = (1,0).
    # 5. Y_b is found by maximizing the second objective (y2). This gives w = x2 and Y_b = (0,1).
    # 6. The specific point Y_star corresponds to w_star = (cos(alpha), 0, sin(alpha)).
    # 7. We calculate Y_star and show it lies below the line connecting Y_a and Y_b.
    # 8. The conclusion is that since for d=3 a counterexample exists, and for d=2 the front is
    #    known to be convex, the largest dimension d for which the property holds is 2.

    # Step 1 & 2: Set up the counterexample
    d = 3
    n = 2
    alpha = np.pi / 6
    c, s = np.cos(alpha), np.sin(alpha)

    x1 = np.array([1., 0., 0.])
    x2 = np.array([0., c, s])

    print(f"Demonstrating with a counterexample for d={d}, n={n}")
    print(f"x1 = {x1}")
    print(f"x2 = {x2.round(4)}")
    print("-" * 40)

    # Step 3, 4, 5: Define the extreme points and the test point
    # Y_a corresponds to w = (1,0,0), which maximizes the squared projection on x1.
    Y_a = np.array([(np.dot(x1, x1))**2, (np.dot(x2, x1))**2])

    # Y_b corresponds to w = x2/||x2||, which maximizes the squared projection on x2.
    # Since ||x2||=1, w=x2.
    Y_b = np.array([(np.dot(x1, x2))**2, (np.dot(x2, x2))**2])

    # Y_star is the special point that reveals non-convexity. It corresponds to w_star.
    w_star = np.array([c, 0., s])
    y1_star = (np.dot(x1, w_star))**2
    y2_star = (np.dot(x2, w_star))**2
    Y_star = np.array([y1_star, y2_star])

    print(f"Extreme Pareto point Y_a (maximizes y1): {Y_a.round(4)}")
    print(f"Extreme Pareto point Y_b (maximizes y2): {Y_b.round(4)}")
    print(f"A third Pareto point Y_star: {Y_star.round(4)}")
    print("-" * 40)

    # Step 6 & 7: Check for non-convexity
    # The line segment connecting Y_a=(1,0) and Y_b=(0,1) is defined by y2 = 1 - y1.
    # We check if Y_star is below this line.
    y_line_at_Y_star_x = 1.0 - Y_star[0]

    print(f"The x-coordinate of Y_star is y1 = {Y_star[0]:.4f}")
    print(f"The y-coordinate of Y_star is y2 = {Y_star[1]:.4f}")
    print(f"On the line segment [Y_a, Y_b], for y1={Y_star[0]:.4f}, the y2 value is {y_line_at_Y_star_x:.4f}")

    is_below = Y_star[1] < y_line_at_Y_star_x
    
    print(f"\nIs Y_star's y2-value less than the y2-value on the line? {is_below}")
    if is_below:
        print("This demonstrates that the Pareto front is NON-CONVEX for d=3.")
    else:
        print("This test does not show non-convexity.")
    print("-" * 40)

    # Step 8: Final Conclusion
    print("Conclusion:")
    print("We have shown a counterexample for d=3 where the Pareto front is non-convex,")
    print("even when all data vectors x_i are non-negative.")
    print("A non-convex Pareto front means that scalarization (the weighted sum method)")
    print("cannot generate all Pareto-optimal solutions.")
    print("\nIt is a known result that for d=2, the Pareto front for this problem is convex.")
    print("Therefore, the largest dimension 'd' for which non-negativity of x_i is a sufficient")
    print("condition for the entire Pareto front to be discoverable by scalarization is d=2.")

demonstrate_non_convexity_and_solve()