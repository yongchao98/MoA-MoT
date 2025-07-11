import numpy as np

def run_simulation():
    """
    This script evaluates the optimal resource distribution strategy for a mother bird
    by simulating different scenarios. It demonstrates that the concavity or convexity
    of the survival function s(r) is the determining factor for the optimal strategy.
    """

    # --- Parameters ---
    n = 5      # Number of offspring
    R = 10.0   # Total resources
    r_max = 4.0  # Max resources per offspring

    # --- Strategies ---
    def fair_strategy(n, R, s_func, s_func_name):
        """Evenly divide all resources among n offspring."""
        r_i = R / n
        total_survival = n * s_func(r_i)
        # Final Equation Output
        print(f"Fair Strategy for {s_func_name}: ", end="")
        for i in range(n):
            print(f"s({r_i:.2f})", end=" + " if i < n - 1 else f" = {total_survival:.4f}\n")
        return total_survival

    def unfair_strategy(n, R, r_max, s_func, s_func_name):
        """Give r_max to k offspring, remainder to one, and zero to the rest."""
        k = int(R / r_max)
        r_rem = R - k * r_max
        zeros = n - k - 1
        
        total_survival = k * s_func(r_max) + s_func(r_rem) + zeros * s_func(0)
        
        # Final Equation Output
        print(f"Unfair Strategy for {s_func_name}: ", end="")
        parts = []
        if k > 0:
            parts.append(f"{k}*s({r_max:.2f})")
        if r_rem > 1e-9: # Only show if remainder is non-trivial
            parts.append(f"1*s({r_rem:.2f})")
        if zeros > 0:
            parts.append(f"{zeros}*s(0.00)")
        
        equation_str = " + ".join(parts)
        print(f"{equation_str} = {total_survival:.4f}\n")
        return total_survival

    # --- Survival Functions s(r) ---
    s_funcs = {
        "Concave Increasing [s(r)=sqrt(r)]": lambda r: np.sqrt(r + 1e-9),
        "Concave Decreasing [s(r)=sqrt(5-r)]": lambda r: np.sqrt(5 - r + 1e-9),
        "Convex Increasing [s(r)=r^2]": lambda r: r**2,
        "Convex Decreasing [s(r)=(r-5)^2]": lambda r: (r - 5)**2
    }
    
    print(f"--- Running Simulations (n={n}, R={R}, r_max={r_max}) ---\n")
    
    results = {}
    # Run simulation for each function type
    for name, func in s_funcs.items():
        print(f"--- Testing with: {name} ---")
        fair_val = fair_strategy(n, R, func, name)
        unfair_val = unfair_strategy(n, R, r_max, func, name)
        results[name] = "Fair" if fair_val > unfair_val else "Unfair"

    print("\n--- Conclusion from Simulations ---")
    # Evaluate statements based on simulation results
    # 1. Test Statement [1]: "increasing -> fair is optimal"
    print(f"1. For 'Concave Increasing', optimal is '{results['Concave Increasing [s(r)=sqrt(r)]']}'.")
    print(f"   For 'Convex Increasing', optimal is '{results['Convex Increasing [s(r)=r^2]']}'.")
    print("   Since an increasing function can have an Unfair optimal strategy, Statement [1] is FALSE.\n")

    # 2. Test Statement [2]: "decreasing -> unfair is optimal"
    print(f"2. For 'Concave Decreasing', optimal is '{results['Concave Decreasing [s(r)=sqrt(5-r)]']}'.")
    print("   Since a decreasing function can have a Fair optimal strategy, Statement [2] is FALSE.\n")

    # 3. Test Statement [3]: "concave increasing -> fair; concave decreasing -> unfair"
    print(f"3. For 'Concave Decreasing', optimal is '{results['Concave Decreasing [s(r)=sqrt(5-r)]']}', not Unfair.")
    print("   Therefore, Statement [3] is FALSE.\n")

    # 4. Test Statement [4]: "concave -> fair is optimal"
    print(f"4. For 'Concave Increasing', optimal is '{results['Concave Increasing [s(r)=sqrt(r)]']}'.")
    print(f"   For 'Concave Decreasing', optimal is '{results['Concave Decreasing [s(r)=sqrt(5-r)]']}'.")
    print("   In all simulated concave cases, Fair is optimal. This supports that Statement [4] is TRUE.\n")

run_simulation()