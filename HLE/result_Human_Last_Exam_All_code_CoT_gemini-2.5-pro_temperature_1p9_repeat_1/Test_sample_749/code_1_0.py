import sympy

def solve_branching_walk():
    """
    Solves the branching random walk problem by analyzing its drift and transience.
    """
    
    # Define the parameter h as a symbolic variable
    h = sympy.Symbol('h')

    # Given jump probabilities
    p_left_red = sympy.Rational(4, 5)
    p_right_red = sympy.Rational(1, 5)
    p_left_blue = sympy.Rational(1, 5)
    p_right_blue = sympy.Rational(4, 5)

    print("Step 1: Formalize the problem.")
    print("The goal is to calculate lim_{h->0} P(site 0 is visited by infinitely many particles).")
    print("Let the event be I_0. We want to find lim_{h->0} P(I_0).")
    print("-" * 50)

    print("Step 2: Analyze the average drift of a single particle.")
    print("The drift depends on the site's color. A site is Red with probability h, and Blue with probability (1-h).")
    
    # Calculate drift (expected displacement) for each color
    # E[jump|Color] = (+1)*p_right + (-1)*p_left
    drift_red = p_right_red - p_left_red
    drift_blue = p_right_blue - p_left_blue
    
    # The overall average drift 'v' is the expectation over the environment
    v = h * drift_red + (1 - h) * drift_blue
    v_simplified = sympy.simplify(v)
    
    print(f"Drift on a Red site: {drift_red}")
    print(f"Drift on a Blue site: {drift_blue}")
    print("By averaging over the probabilities of a site being Red or Blue, we find the mean drift v(h):")
    print(f"v(h) = h * ({drift_red}) + (1-h) * ({drift_blue}) = {v_simplified}")
    
    # Analyze the drift in the limit h -> 0
    v_limit = sympy.limit(v_simplified, h, 0)
    print(f"\nAs h approaches 0, the drift approaches: lim_{h->0} v(h) = {v_limit}")
    print("Since the limiting drift is positive, for any small h>0, the particle system has a strong net drift to the right.")
    print("-" * 50)
    
    print("Step 3: Confirm transience using the formal criterion for RWRE.")
    print("A Random Walk in a Random Environment (RWRE) on Z is transient if E[log(rho)] is non-zero,")
    print("where rho is the ratio of left-to-right jump probabilities at a site.")
    
    rho_red = p_left_red / p_right_red
    rho_blue = p_left_blue / p_right_blue
    print(f"\nrho_red = P(Left|Red) / P(Right|Red) = {p_left_red} / {p_right_red} = {rho_red}")
    print(f"rho_blue = P(Left|Blue) / P(Right|Blue) = {p_left_blue} / {p_right_blue} = {rho_blue}")
    
    # The expectation of log(rho)
    E_log_rho = h * sympy.log(rho_red) + (1 - h) * sympy.log(rho_blue)
    E_log_rho_simplified = sympy.simplify(E_log_rho)
    E_log_rho_limit = sympy.limit(E_log_rho_simplified, h, 0)
    
    print("\nThe expected value of log(rho) is E[log(rho)](h):")
    # Using unicode for rho for better display
    print(f"E[log(\u03C1)](h) = h*log({rho_red}) + (1-h)*log({rho_blue}) = {E_log_rho_simplified}")
    print(f"\nAs h approaches 0, this value approaches: lim_{h->0} E[log(\u03C1)](h) = {E_log_rho_limit.evalf()} ({E_log_rho_limit})")

    print("\nSince E[log(rho)] is non-zero for h in (0, 1/2), the underlying walk is transient.")
    print("A negative value indicates a strong drift to the right (as p_right > p_left).")
    print("-" * 50)

    print("Step 4: Relate the single-particle walk to the entire branching process.")
    print("A key theorem in the study of Branching Random Walks states that the branching process is transient")
    print("if the underlying random walk is transient. This means the entire cloud of particles will move together")
    print("in the direction of the drift.")
    print("\nOur process starts at site 3. With a strong drift to the right (towards +infinity), the entire particle")
    print("population will almost surely move away from the origin.")
    print("-" * 50)

    print("Step 5: Conclusion and Final Calculation.")
    print("Because the particle cloud drifts to +infinity, any fixed site (like site 0) will eventually")
    print("be left behind by the entire population. Therefore, site 0 can only be visited by a finite number of particles.")
    print("\nThis means that the probability of site 0 being visited by infinitely many particles is 0 for any h in (0, 1/2).")
    print("P(I_0) = 0 for h \u2208 (0, 1/2)")

    final_result = 0
    print("\nThe limit of a function that is constantly 0 is 0. So the final equation and result is:")
    print(f"lim_{h->0} P(I_0) = {final_result}")

if __name__ == "__main__":
    solve_branching_walk()
    print("\nFinal Answer")
    print(f"The equation is lim_{h->0} P(I_0) = 0. The numbers in this equation are 0.")