def solve_cohomology_rank():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology ring
    of the space of 3-subspaces of R^5.
    """

    # Define the parameters of the Grassmannian
    k = 3
    n = 5

    # Step 1: Explain the space and the object of study.
    print(f"The problem asks for the rank of the torsion subgroup of the integral cohomology ring of the space of {k}-subspaces of R^{n}.")
    print(f"This space is the real Grassmannian manifold Gr({k}, {n}).")
    print(f"We are interested in the group G = H*(Gr({k}, {n}); Z).")
    print("-" * 20)

    # Step 2: State the key mathematical fact.
    print("A fundamental theorem in algebraic topology states that the integral cohomology ring of a real Grassmannian, H*(Gr(k, n); Z), is a free abelian group.")
    print("This is because the Grassmannian has a cell decomposition (the Schubert cells) which leads to a cellular chain complex with free homology groups.")
    print("-" * 20)

    # Step 3: Deduce the structure of the torsion subgroup.
    print("A free abelian group, by definition, has no non-zero elements of finite order.")
    print("The torsion subgroup of a group consists of all its elements of finite order.")
    print("Therefore, the torsion subgroup of H*(Gr(3, 5); Z) is the trivial group, {0}.")
    print("-" * 20)

    # Step 4: Calculate the rank.
    print("The rank of an abelian group is the number of Z summands in its free part.")
    print("The torsion subgroup is {0}, which is a finite group. A finite group has no free part.")
    final_rank = 0
    print(f"Thus, its rank is {final_rank}.")
    print("-" * 20)

    # Step 5: Display the final equation as requested.
    print("Final Equation:")
    print(f"rank(TorsionSubgroup(H*(Gr({k}, {n}); Z))) = {final_rank}")

if __name__ == "__main__":
    solve_cohomology_rank()