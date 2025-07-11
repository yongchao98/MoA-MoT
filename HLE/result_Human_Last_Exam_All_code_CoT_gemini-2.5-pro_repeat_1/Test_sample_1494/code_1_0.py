import numpy as np
import scipy.optimize

def analyze_decay_scaling(t_values, a, c):
    """
    Analyzes the scaling of the dominant wavenumber and the decay exponent.
    
    Args:
        t_values (list or np.array): A list of time points to analyze.
        a (float): Parameter for the initial data spectrum decay (e.g., S_k ~ exp(-a*k)).
        c (float): Parameter for the modal decay rate (mu_k ~ c/k^2).
    """
    
    print(f"Analysis with parameters a={a}, c={c}\n")
    print(f"{'Time (t)':<12} | {'Optimal k*':<15} | {'k* / t^(1/3)':<15} | {'Exponent f(k*)':<20} | {'f(k*) / t^(1/3)':<15}")
    print("-" * 85)

    k_stars = []
    f_mins = []

    for t in t_values:
        # The function to minimize is the exponent: f(k) = a*k + 2*c*t/k^2
        def exponent_func(k):
            # We use a large penalty for k<=0 to ensure k is positive
            if k <= 0:
                return 1e20
            return a * k + 2 * c * t / (k**2)

        # Find the k that minimizes the exponent function
        # We start the search from a reasonable guess, e.g., t^(1/3)
        result = scipy.optimize.minimize(exponent_func, x0=(t+1)**(1/3), bounds=[(1e-9, None)])
        k_star = result.x[0]
        f_min = result.fun
        
        k_stars.append(k_star)
        f_mins.append(f_min)

        # Theoretical prediction for the ratio k* / t^(1/3)
        k_ratio = k_star / (t**(1/3)) if t > 0 else 0
        
        # Theoretical prediction for the ratio f(k*) / t^(1/3)
        f_ratio = f_min / (t**(1/3)) if t > 0 else 0
        
        print(f"{t:<12.1e} | {k_star:<15.4f} | {k_ratio:<15.4f} | {f_min:<20.4f} | {f_ratio:<15.4f}")

    # You can also use these lists for plotting if in a notebook
    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(12, 5))
    # plt.subplot(1, 2, 1)
    # plt.loglog(t_values, k_stars, 'o-', label='Numerical k*')
    # plt.loglog(t_values, [t**(1/3) for t in t_values], '--', label='t^(1/3) scaling')
    # plt.title('Optimal Wavenumber k* vs. Time')
    # plt.xlabel('Time (t)')
    # plt.ylabel('k*')
    # plt.legend()
    # plt.grid(True)
    # plt.subplot(1, 2, 2)
    # plt.loglog(t_values, f_mins, 'o-', label='Numerical min exponent')
    # plt.loglog(t_values, [t**(1/3) for t in t_values], '--', label='t^(1/3) scaling')
    # plt.title('Minimal Exponent f(k*) vs. Time')
    # plt.xlabel('Time (t)')
    # plt.ylabel('f(k*)')
    # plt.legend()
    # plt.grid(True)
    # plt.tight_layout()
    # plt.show()


# --- Main execution ---
if __name__ == "__main__":
    # Parameters for the analysis
    a_param = 1.0  # Decay of initial spectrum
    c_param = 1.0  # Strength of high-k damping

    # A range of time values, spanning several orders of magnitude
    time_points = np.logspace(2, 8, num=7)
    
    analyze_decay_scaling(time_points, a_param, c_param)
