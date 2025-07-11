import sympy

def display_proof():
    """
    This function prints a step-by-step proof determining if the solution
    to the given Cauchy problem can blow-up in finite time.
    """

    # --- Introduction and Problem Statement ---
    print("The question is whether the solution to the following Cauchy problem on R^3 can blow up in finite time from smooth divergence-free initial data u_0.")
    print("\nProblem:")
    print("  ∂_t u + u·∇u + (1+t)Δu - ∇p = 0")
    print("  ∇·u = 0")
    print("  u(x, t=0) = u_0(x)")
    print("-" * 50)
    print("\nAnswer: No, the solution cannot blow up in finite time.\n")
    print("Here is the proof by contradiction:")
    print("-" * 50)

    # --- Step 1: L^2 Energy Estimate ---
    print("\nStep 1: Derive the L^2 energy estimate.")
    print("We take the L^2 inner product of the momentum equation with u. This means we multiply the equation by u and integrate over all of space (R^3).")
    print("The resulting equation is:")
    print("  ∫(∂_t u)·u dx + ∫(u·∇u)·u dx + ∫(1+t)Δu·u dx - ∫∇p·u dx = 0\n")

    print("Let's analyze each term:")
    print("1. ∫(∂_t u)·u dx = (1/2) d/dt ∫|u|^2 dx = (1/2) d/dt ||u||_L2^2")
    print("   This term represents the rate of change of the kinetic energy of the fluid.")
    print("2. ∫(u·∇u)·u dx = 0")
    print("   The nonlinear term vanishes for divergence-free fields (∇·u = 0), assuming u decays at infinity.")
    print("3. -∫∇p·u dx = 0")
    print("   The pressure term vanishes after integration by parts, again using ∇·u = 0.")
    print("4. ∫(1+t)Δu·u dx = -(1+t) ∫|∇u|^2 dx = -(1+t) ||∇u||_L2^2")
    print("   This term, after integration by parts, represents the energy dissipation due to viscosity.\n")

    # --- Step 2: The Energy Equality ---
    print("Step 2: State the resulting energy equality.")
    print("Combining the terms, we get a differential equation for the kinetic energy:")
    print("  (1/2) d/dt ||u(t)||_L2^2 + (1+t) ||∇u(t)||_L2^2 = 0")
    print("This can be rewritten as:")
    print("  d/dt ||u(t)||_L2^2 = -2(1+t) ||∇u(t)||_L2^2\n")

    # --- Step 3: Integrate the Energy Equality ---
    print("Step 3: Integrate the energy equality in time.")
    print("We integrate the equation from t=0 to an arbitrary time T > 0:")
    print("  ∫[0, T] d/dt ||u(t)||_L2^2 dt = -2 ∫[0, T] (1+t) ||∇u(t)||_L2^2 dt")
    print("  ||u(T)||_L2^2 - ||u_0||_L2^2 = -2 ∫[0, T] (1+t) ||∇u(t)||_L2^2 dt\n")
    print("Rearranging and using the fact that ||u(T)||_L2^2 ≥ 0, we obtain a crucial bound:")
    print("  ∫[0, T] (1+t) ||∇u(t)||_L2^2 dt = (1/2) (||u_0||_L2^2 - ||u(T)||_L2^2) ≤ (1/2) ||u_0||_L2^2\n")
    print("This inequality holds for any time T > 0, which means the total weighted dissipated energy is bounded by the initial kinetic energy.")

    # --- Step 4: Proof by Contradiction ---
    print("\nStep 4: The Proof by Contradiction.")
    print("Now, let's assume for the sake of contradiction that the solution *does* blow up at a finite time T_c > 0.")
    print("By definition, a 'blow-up' for a smooth solution means that the H^1 norm (or a higher norm) becomes infinite at T_c.")
    print("Specifically, the L^2 norm of the gradient becomes unbounded:")
    print("  lim_{t → T_c⁻} ||∇u(t)||_L2^2 = ∞\n")

    print("Let's consider the integral we bounded in Step 3, up to the blow-up time T_c:")
    print("  I = ∫[0, T_c] (1+t) ||∇u(t)||_L2^2 dt\n")
    print("The integrand, f(t) = (1+t) ||∇u(t)||_L2^2, is a continuous and positive function on the interval [0, T_c).")
    print("As t approaches T_c from below:")
    print("  - The term (1+t) approaches a finite positive value (1+T_c).")
    print("  - The term ||∇u(t)||_L2^2 approaches infinity (by our assumption).")
    print("Therefore, the integrand f(t) must also approach infinity as t → T_c⁻.\n")
    print("A positive function that approaches infinity at the endpoint of a finite interval must have an infinite integral over that interval.")
    print("So, based on our assumption of a finite-time blow-up, we must have:")
    print("  I = ∫[0, T_c] (1+t) ||∇u(t)||_L2^2 dt = ∞\n")

    # --- Step 5: The Contradiction ---
    print("Step 5: Reaching the Contradiction.")
    print("Let's summarize our findings:")
    print("1. From the energy equality (Step 3), we proved: ∫[0, T_c] (1+t) ||∇u(t)||_L2^2 dt ≤ (1/2) ||u_0||_L2^2")
    print("2. From the assumption of a finite-time blow-up (Step 4), we deduced: ∫[0, T_c] (1+t) ||∇u(t)||_L2^2 dt = ∞")

    # Use a hypothetical number for the initial energy for a clear final statement.
    initial_energy_sq = 10.0
    bound = initial_energy_sq / 2.0
    print(f"\nLet's assume the initial energy ||u_0||_L2^2 is {initial_energy_sq}. The bound from Step 3 becomes {bound}.")

    print("\nThis leads to the following contradictory statement:")
    inf = sympy.oo
    # This print statement fulfills the "output each number in the final equation" requirement.
    print(f'  {inf} <= {bound}')

    print("\nThis is a logical impossibility. Therefore, our initial assumption—that the solution can blow up in finite time—must be false.")

if __name__ == '__main__':
    display_proof()
    print("\n<<<No>>>")