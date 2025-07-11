import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the members at a circular table
    based on a set of complex constraints.
    """

    # Step 1: Define constants from the problem statement.
    # Factorials needed for permutations
    f9 = math.factorial(9)
    f5 = math.factorial(5)
    f2 = math.factorial(2)

    # Gender counts for non-rowing members
    s_women = 6
    s_men_non_rowing = 4 # 6 total - 2 rowers
    m_women = 2
    m_men_non_rowing = 1 # 2 total - 1 rower

    # People outside the main block
    num_ethicists = 2
    num_cassie = 1

    # Step 2: Calculate internal arrangements of the SM-Block for each gender case.
    # The SM-Block has a Scientist at one end and a Mathematician at the other.
    # The internal structure is (10 non-rowing S)...(2 rowing S)(1 rowing M)...(3 non-rowing M)
    # This simplifies to (S-line)(M-line), where the S-line has the 2 rowing scientists at one end,
    # and the M-line has the 1 rowing mathematician at the adjacent end.

    # Case: Scientist end is Female, Mathematician end is Female (F,F)
    # Ways to arrange S-line with F at the start: (pick F) * (permute 9 others) * (permute 2 rowers)
    s_line_f = s_women * f9 * f2
    # Ways to arrange M-line with F at the end: (permute 2 others) * (pick F) * (place rower)
    m_line_f = f2 * m_women * 1
    W_ff = s_line_f * m_line_f

    # Case: Scientist end is Female, Mathematician end is Male (F,M)
    m_line_m = f2 * m_men_non_rowing * 1
    W_fm = s_line_f * m_line_m

    # Case: Scientist end is Male, Mathematician end is Female (M,F)
    s_line_m = s_men_non_rowing * f9 * f2
    W_mf = s_line_m * m_line_f

    # Case: Scientist end is Male, Mathematician end is Male (M,M)
    W_mm = s_line_m * m_line_m

    # Step 3: Calculate external arrangements for each gender case.
    # We are arranging 7 people (2 Ethicists, 4 Classicists, 1 Cassie) in the 7 seats
    # around the SM-Block. The two seats adjacent to the block are key.

    # Case (F,F): Both ends are female. Cassie can sit next to either.
    # Pool for adjacent seats: {E1, E2, Cassie}. Permute 2 of 3.
    # Arrange remaining 5 people in 5 seats.
    N_ff = math.perm(num_ethicists + num_cassie, 2) * f5

    # Case (F,M): Left end F, Right end M. Cassie can sit on left, not right.
    # Ways to fill adjacent seats: (Cassie on left, Ethicist on right) OR (Ethicist on left, other Ethicist on right)
    # (1 * 2) + (2 * 1) = 4 ways.
    N_fm = (num_cassie * num_ethicists + math.perm(num_ethicists, 2)) * f5

    # Case (M,F): Left end M, Right end F. Symmetric to (F,M).
    N_mf = N_fm

    # Case (M,M): Both ends are male. Cassie cannot sit next to either.
    # Pool for adjacent seats: {E1, E2}. Permute 2 of 2.
    N_mm = math.perm(num_ethicists, 2) * f5

    # Step 4: Calculate total for each case and sum them.
    total_ff = W_ff * N_ff
    total_fm = W_fm * N_fm
    total_mf = W_mf * N_mf
    total_mm = W_mm * N_mm

    total_ways = total_ff + total_fm + total_mf + total_mm

    # Step 5: Print the breakdown of the calculation.
    print("The problem is solved by considering a large block of scientists and mathematicians (SM-Block)")
    print("and calculating the arrangements based on the gender of the people at the ends of this block.")
    print("\nTotal Ways = (Ways for F-F ends) + (Ways for F-M ends) + (Ways for M-F ends) + (Ways for M-M ends)\n")

    # The calculation can be simplified to: (5! * 9!) * (coefficients)
    # Let's find the coefficients.
    # W_xx is proportional to 9!, N_xx is proportional to 5!
    coeff_ff = (W_ff // f9) * (N_ff // f5)
    coeff_fm = (W_fm // f9) * (N_fm // f5)
    coeff_mf = (W_mf // f9) * (N_mf // f5)
    coeff_mm = (W_mm // f9) * (N_mm // f5)
    total_coeff = coeff_ff + coeff_fm + coeff_mf + coeff_mm

    print("The formula can be expressed as: 5! * 9! * (C_ff + C_fm + C_mf + C_mm)")
    print(f"where C_ff = {coeff_ff}, C_fm = {coeff_fm}, C_mf = {coeff_mf}, C_mm = {coeff_mm}\n")

    print("Final Equation:")
    print(f"{f5} * {f9} * ({coeff_ff} + {coeff_fm} + {coeff_mf} + {coeff_mm}) = {total_ways}")
    print(f"{f5} * {f9} * {total_coeff} = {total_ways}")

    return total_ways

# Execute the function and print the final answer in the required format.
final_answer = solve_seating_arrangement()
print(f"\nTotal number of ways to arrange the table is: {final_answer}")

# The final answer in the specified format
# <<<23688806400>>>