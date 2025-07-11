def solve_cost_calculation():
    """
    Calculates the minimum cost for computing 2A - 3B on a twisted Edwards curve.

    The plan is to find an optimal sequence of operations and sum their costs,
    based on standard formulas from the Explicit-Formulas Database (EFD),
    assuming the cost of a squaring is the same as a multiplication (S=M).

    The chosen computation path is: 2A - 3B  =>  2(A - B) - B  =>  2(A + (-B)) + (-B)
    """

    # Define costs of primitive operations from the EFD for twisted Edwards curves
    # in extended coordinates (https://www.hyperelliptic.org/EFD/g1p/auto-twisted-extended.html)
    # assuming S=M.

    # Cost to convert an affine point (x,y) to extended (X,Y,Z,T) is 1M for T=x*y.
    cost_affine_to_extended = 1

    # Cost for mixed-coordinate addition of an extended and an affine point (add-2008-hwcd-3).
    # This is unified and costs 8M.
    cost_mixed_add = 8

    # Cost for "double-and-add" of two extended points (dadd-2008-hwcd).
    # Cost is 6M + 4S. With S=M, this is 6M + 4M = 10M.
    cost_dadd_extended = 10

    # --- Calculation for 2A - 3B, rewritten as 2(A + (-B)) + (-B) ---

    # Stage 1: Compute P = A + (-B)
    # Negating B is free. To add affine A and affine -B to get extended P:
    # 1a. Convert A (affine) to A_ext (extended).
    cost_step_1a = cost_affine_to_extended
    # 1b. Add A_ext and -B (affine) using mixed addition.
    cost_step_1b = cost_mixed_add
    cost_stage_1 = cost_step_1a + cost_step_1b

    # Stage 2: Compute R = 2P + (-B)
    # We have P (extended) and -B (affine). To use the efficient extended dadd formula:
    # 2a. Convert -B (affine) to -B_ext (extended).
    cost_step_2a = cost_affine_to_extended
    # 2b. Compute 2P + (-B)_ext using the extended double-and-add formula.
    cost_step_2b = cost_dadd_extended
    cost_stage_2 = cost_step_2a + cost_step_2b

    # Total cost is the sum of all steps.
    total_cost = cost_stage_1 + cost_stage_2

    # Print the detailed breakdown and the final result.
    print("To compute 2A - 3B with the smallest cost, we rewrite the operation as 2(A - B) - B.")
    print("The cost is calculated in two main stages:\n")

    print(f"Stage 1: Compute P = A - B")
    print(f"  - Convert point A (affine) to extended coordinates. Cost: {cost_step_1a}M")
    print(f"  - Perform a mixed addition of A_ext and -B (affine). Cost: {cost_step_1b}M")
    print(f"  - Subtotal for Stage 1: {cost_stage_1}M\n")

    print(f"Stage 2: Compute the final result R = 2P - B")
    print(f"  - Convert point -B (affine) to extended coordinates. Cost: {cost_step_2a}M")
    print(f"  - Perform a 'double-and-add' operation 2P + (-B)_ext. Cost: {cost_dadd_extended}M")
    print(f"  - Subtotal for Stage 2: {cost_stage_2}M\n")

    print("Total Cost Calculation:")
    print(f"Total Cost = (Cost of Stage 1) + (Cost of Stage 2)")
    print(f"Total Cost = ({cost_step_1a}M + {cost_step_1b}M) + ({cost_step_2a}M + {cost_dadd_extended}M)")
    print(f"Total Cost = {cost_stage_1}M + {cost_stage_2}M = {total_cost}M")
    print("\nThe smallest cost is 20 multiplications.")

solve_cost_calculation()