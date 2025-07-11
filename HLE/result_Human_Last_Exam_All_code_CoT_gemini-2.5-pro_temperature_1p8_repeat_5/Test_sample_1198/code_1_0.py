import decimal

def demonstrate_search(mod, n_max=30, a_range=[1, 5]):
    """
    This function numerically demonstrates the existence (or non-existence)
    of a number 'a' for a given modulus by iteratively refining the possible
    intervals for 'a'.
    """
    
    # Use high-precision decimals to avoid floating point issues
    decimal.getcontext().prec = 50

    def n_root(d, n):
        """Calculates the nth root of a Decimal d."""
        if d < 0: return -((-d).ln() / decimal.Decimal(n)).exp()
        if d == 0: return decimal.Decimal(0)
        return (d.ln() / decimal.Decimal(n)).exp()

    def get_good_intervals(current_intervals, n, mod):
        """Filters the current intervals to find sub-intervals satisfying the condition for n."""
        next_intervals = []
        required_rem = n % mod

        for min_a, max_a in current_intervals:
            min_an = min_a ** n
            max_an = max_a ** n

            start_k = min_an.to_integral_value(rounding=decimal.ROUND_FLOOR)
            end_k = max_an.to_integral_value(rounding=decimal.ROUND_CEILING)

            k = start_k
            while k < end_k:
                if k % mod == required_rem:
                    # Condition: k <= a^n < k+1
                    # Interval for 'a': k^(1/n) <= a < (k+1)^(1/n)
                    k_dec = decimal.Decimal(k)
                    new_min = n_root(k_dec, n)
                    new_max = n_root(k_dec + 1, n)

                    # Intersect with the parent interval
                    overlap_min = max(new_min, min_a)
                    overlap_max = min(new_max, max_a)

                    if overlap_min < overlap_max:
                        next_intervals.append([overlap_min, overlap_max])
                k += 1
        return next_intervals

    print(f"\n--- Demonstrating search for modulus m = {mod} ---")
    
    candidate_intervals = [[decimal.Decimal(x) for x in a_range]]

    for n in range(1, n_max + 1):
        if not candidate_intervals:
            print(f"n={n}: No candidate intervals left.")
            break
        
        candidate_intervals = get_good_intervals(candidate_intervals, n, mod)
        
        # Output summary to avoid printing too much data
        print(f"n={n}: Found {len(candidate_intervals)} candidate interval(s).")
        if candidate_intervals:
            # Show a representative sample of the total range
            total_min = candidate_intervals[0][0]
            total_max = candidate_intervals[-1][1]
            print(f"       Total range of candidates spans approx: [{total_min:.6f}, {total_max:.6f})")


    print("-" * 30)
    if candidate_intervals:
        print(f"Conclusion for mod={mod}: A solution likely exists.")
    else:
        print(f"Conclusion for mod={mod}: A solution does NOT seem to exist.")
    print("-" * 30)


# --- Run Demonstrations ---
# For modulo 2, the number of valid intervals remains stable.
demonstrate_search(mod=2)

# For modulo 3, the search fails as the set of candidates becomes empty.
demonstrate_search(mod=3, a_range=[1,2])

print("\nFinal Answer: Yes, No")
