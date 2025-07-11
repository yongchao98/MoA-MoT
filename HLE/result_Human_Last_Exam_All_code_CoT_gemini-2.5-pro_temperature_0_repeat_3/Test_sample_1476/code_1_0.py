import numpy as np
from scipy.linalg import null_space

def solve():
    """
    This function explains the logical deduction and provides a coded example.
    """
    print("### Step-by-step Deduction ###")
    print("1. The condition 'no cycles having non-zero sum' means the signal x¹ is orthogonal to the cycle space, ker(B₁).")
    print("   This implies x¹ must be in the orthogonal complement of ker(B₁), which is im(B₁ᵀ), the space of gradients.")
    print("\n2. The condition 'B₁x¹1ᵀ = 0' implies that B₁x¹ = 0.")
    print("   This means x¹ is in the kernel of B₁, i.e., x¹ ∈ ker(B₁).")
    print("\n3. So, x¹ is in both a subspace (im(B₁ᵀ)) and its orthogonal complement (ker(B₁)).")
    print("   The only vector that satisfies this is the zero vector. Therefore, we must conclude that x¹ = 0.")
    print("\n4. The total variation (TV) of the vertex signal x⁰ is defined as TV = Σ |x⁰ᵤ - x⁰ᵥ| = Σ x¹ₑ.")
    print("   Since x¹ = 0, the total variation must be 0.")
    print("   This corresponds to answer choice D.")

    print("\n### Code Demonstration ###")
    print("Let's use a triangle graph as an example, which has one cycle.")
    # V = {0, 1, 2}, Edges (oriented): e0=(0,1), e1=(1,2), e2=(2,0)
    B1 = np.array([
        [-1, 0, 1],
        [1, -1, 0],
        [0, 1, -1]
    ])
    B1_T = B1.T
    
    # The basis for the cycle space ker(B₁) for a triangle is a vector representing the cycle.
    # For the orientation 0->1->2->0, the cycle vector is [1, 1, 1]
    ker_B1_basis = null_space(B1)
    # Ensure the basis vector is simple for demonstration
    if ker_B1_basis.shape[1] > 0:
        c = ker_B1_basis[:, 0] * np.sign(ker_B1_basis[0,0])
    else:
        c = None
    print(f"\nIncidence Matrix B₁:\n{B1}")
    if c is not None:
        print(f"\nA basis vector for the cycle space ker(B₁) is c = {c.flatten().round(2)}")

    print("\n--- Case 1: A signal that meets the problem's conditions ---")
    # If x¹=0, then |x⁰ᵤ - x⁰ᵥ|=0 for all edges, so x⁰ must be constant on connected components.
    x0_const = np.array([5, 5, 5])
    print(f"Let's use a constant vertex signal x⁰ = {x0_const}")
    
    grad_x0 = B1_T @ x0_const
    print(f"Gradient of x⁰ is B₁ᵀx⁰ = {grad_x0}")
    x1 = np.abs(grad_x0)
    print(f"Resulting edge signal x¹ = |B₁ᵀx⁰| = {x1}")

    # Check Premise 1: x¹ is orthogonal to ker(B₁)
    if c is not None:
        dot_product = x1.T @ c
        print(f"Checking Premise 1 (orthogonality): cᵀx¹ = {c.flatten().round(2)} . {x1} = {dot_product:.2f}")
        print(f"This is 0, so Premise 1 holds.")

    # Check Premise 2: x¹ is in ker(B₁)
    div_x1 = B1 @ x1
    print(f"Checking Premise 2 (divergence-free): B₁x¹ = \n{B1}\n @ {x1} = {div_x1}")
    print(f"This is the zero vector, so Premise 2 holds.")

    # Calculate Total Variation
    total_variation = np.sum(x1)
    print(f"\nFor this case, both premises hold. The total variation is Σ x¹ₑ = {total_variation}.")
    print("This confirms our deduction: when the premises hold, the total variation is 0.")

    print("\n--- Case 2: A signal that does NOT meet the conditions ---")
    x0_nonconst = np.array([1, 2, 5])
    print(f"Let's use a non-constant vertex signal x⁰ = {x0_nonconst}")
    
    grad_x0 = B1_T @ x0_nonconst
    print(f"Gradient of x⁰ is B₁ᵀx⁰ = {grad_x0}")
    x1 = np.abs(grad_x0)
    print(f"Resulting edge signal x¹ = |B₁ᵀx⁰| = {x1}")

    # Check Premise 1: x¹ is orthogonal to ker(B₁)
    if c is not None:
        dot_product = x1.T @ c
        print(f"Checking Premise 1 (orthogonality): cᵀx¹ = {c.flatten().round(2)} . {x1} = {dot_product:.2f}")
        print(f"This is not 0, so Premise 1 fails.")

    # Check Premise 2: x¹ is in ker(B₁)
    div_x1 = B1 @ x1
    print(f"Checking Premise 2 (divergence-free): B₁x¹ = \n{B1}\n @ {x1} = {div_x1}")
    print(f"This is not the zero vector, so Premise 2 also fails.")
    
    total_variation = np.sum(x1)
    print(f"\nIn this case, the premises do not hold, and the total variation is {total_variation}, which is not 0.")

solve()
print("<<<D>>>")