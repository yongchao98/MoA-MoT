def solve_hodge_number():
    """
    Calculates the maximal value of the Hodge number h^{1,1}(M).
    
    M is the smooth manifold obtained by blowing up the singular locus of (S x C) / (rho x psi),
    where S is a K3 surface, C is a genus 2 curve, rho is a non-symplectic involution on S,
    and psi is an involution on C.
    """

    # Step 1: Calculate h^1,1 of the quotient space X = (S x C) / (rho x psi).
    # h^1,1(X) is the dimension of the (+1)-eigenspace of the (rho x psi)^* action
    # on H^1,1(S x C).

    # For a non-symplectic involution rho on a K3 surface S, the dimension of the
    # invariant part of H^2(S) is 10. This part is of type (1,1).
    h11_plus_S = 10
    
    # The involution psi acts trivially on H^1,1(C) and H^0,0(C).
    # h^0,0(S) is also invariant.
    h11_plus_C = 1
    h00_plus_S = 1
    h00_plus_C = 1

    # Using the Künneth formula for the invariant part of the cohomology:
    # h^1,1(X) = h^1,1_+(S x C) = h^1,1_+(S)h^0,0_+(C) + h^0,0_+(S)h^1,1_+(C)
    h11_X = h11_plus_S * h00_plus_C + h00_plus_S * h11_plus_C
    
    print(f"The Hodge number h^1,1 of the quotient space X is: {h11_X}")

    # Step 2: Determine the contribution from the blow-up.
    # The contribution is the number of connected components of the singular locus.
    # Singular locus = Fix(rho) x Fix(psi).
    # Number of components = (components of Fix(rho)) * (components of Fix(psi)).

    # To maximize this, we need to maximize the number of components for each fixed locus.
    
    # For an involution psi on a genus 2 curve C, the maximum number of fixed points
    # is 6 (for the hyperelliptic involution).
    N_psi_max = 6
    print(f"The maximum number of fixed points for psi on C (genus 2) is N_psi = {N_psi_max}.")

    # For a non-symplectic involution rho on a K3 surface S, the maximum number of
    # connected components k of the fixed locus is 10.
    k_max = 10
    print(f"The maximum number of connected components of the fixed locus of rho on S is k = {k_max}.")

    # The number of new classes in H^1,1 from the blow-up is k * N_psi.
    blow_up_contribution = k_max * N_psi_max
    print(f"The contribution to h^1,1 from the blow-up is k * N_psi = {k_max} * {N_psi_max} = {blow_up_contribution}.")

    # Step 3: Calculate the maximal h^1,1(M).
    # h^1,1(M) = h^1,1(X) + (blow-up contribution)
    h11_M_max = h11_X + blow_up_contribution

    # Final result presentation
    print("\nThe maximal value of h^1,1(M) is calculated as:")
    print(f"h^1,1(M) = h^1,1(X) + k_max * N_psi_max")
    # The user requested to output each number in the final equation.
    print(f"h^1,1(M) = {h11_X} + {k_max} * {N_psi_max} = {h11_M_max}")
    
    return h11_M_max

if __name__ == '__main__':
    max_hodge_number = solve_hodge_number()
    # The final answer is returned in the specified format.
    # print(f"\n<<< {max_hodge_number} >>>")