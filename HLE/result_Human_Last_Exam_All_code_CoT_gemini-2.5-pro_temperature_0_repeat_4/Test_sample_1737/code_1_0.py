import sys

def count_d_values(N):
    """
    Calculates the number of different numerical values the non-zero d_ijk for SU(N) take.

    Args:
        N (int): The dimension of the special unitary group SU(N).

    Returns:
        None. Prints the result and the calculation breakdown.
    """
    if not isinstance(N, int) or N < 2:
        print("Please provide an integer N >= 2.", file=sys.stderr)
        return

    # For SU(2), all d_ijk are zero.
    if N == 2:
        print("For SU(2), all d_ijk coefficients are zero.")
        print("Number of non-zero values = 0")
        print("<<<0>>>")
        return

    # For N >= 3, we can calculate the number of values by summing up
    # the new values introduced by each of the N-1 diagonal generators.
    
    # l=1 diagonal generator (and the purely off-diagonal case)
    # contributes 2 values (+-1/2) if N>=3.
    terms = [2]
    
    # l=2 diagonal generator
    if N == 3:
        # For N=3, l=2 is the last generator. It doesn't have access to
        # indices k > l+1, which restricts the number of values.
        terms.append(3)
    elif N > 3:
        # For N>=4, the l=2 generator contributes 4 distinct values.
        terms.append(4)

    # l=3 to N-2 diagonal generators
    # Each contributes 5 new distinct values.
    if N > 4:
        for l in range(3, N - 1):
            terms.append(5)
            
    # l=N-1 (last) diagonal generator (for N>=4)
    # This generator is also restricted and contributes 3 new values.
    if N > 3:
        terms.append(3)

    total_values = sum(terms)
    
    # Building the explanation string
    explanation = f"For SU({N}), the number of distinct non-zero d_ijk values is calculated as follows:\n"
    calculation_str = " + ".join(map(str, terms))
    explanation += f"Number of values = {calculation_str} = {total_values}"
    
    print(explanation)
    print(f"<<<{total_values}>>>")


if __name__ == '__main__':
    try:
        # You can change the value of N here to test for different groups.
        N_input = 5
        print(f"Calculating for N = {N_input}\n")
        count_d_values(N_input)
    except (ValueError, IndexError):
        print("Please provide a valid integer N as an argument.", file=sys.stderr)
        print("Example: python your_script_name.py 5", file=sys.stderr)
