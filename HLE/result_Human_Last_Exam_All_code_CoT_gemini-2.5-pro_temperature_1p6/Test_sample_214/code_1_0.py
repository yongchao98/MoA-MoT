import collections

class PosetCounter:
    """
    A class to perform the necessary computations for the problem.
    It assumes S is the set of divisors of a number n.
    """
    def __init__(self, n):
        self.n = n
        self.divisors_of_n = self._get_divisors(n, all_divs=True)
        self.poset_elements = sorted(list(self.divisors_of_n))

        # Step 1: Compute τ, the set of all upper sets of the divisibility poset.
        self.tau = self._get_upper_sets()
        self.tau_poset_elements = sorted(list(self.tau), key=len)
        
        # Step 2: Count the lower sets of the poset (τ, ⊆).
        self.num_open_sets = self._count_lower_sets_of_tau()

    def _get_divisors(self, num, all_divs=False):
        """Helper to get divisors of a number."""
        divs = set()
        if all_divs:
            limit = int(num**0.5) + 1
            for i in range(1, limit):
                if num % i == 0:
                    divs.add(i)
                    divs.add(num // i)
        else: # Proper divisors within the main set S
            for d in self.poset_elements:
                if d < num and num % d == 0:
                    divs.add(d)
        return frozenset(divs)

    def _get_upper_sets(self):
        """
        Computes the set of all upper sets of the divisibility poset on S.
        This is done by finding all lower sets (divisor-closed sets) and taking complements.
        """
        lower_sets = {frozenset()}
        for k in self.poset_elements:
            new_sets_for_k = set()
            proper_divs = self._get_divisors(k)
            for ls in lower_sets:
                if proper_divs.issubset(ls):
                    new_sets_for_k.add(ls.union({k}))
            lower_sets.update(new_sets_for_k)
        
        full_set = frozenset(self.poset_elements)
        upper_sets = {full_set - ls for ls in lower_sets}
        return upper_sets

    def _count_lower_sets_of_tau(self):
        """
        Counts the number of lower sets of the poset (τ, ⊆).
        This is equivalent to counting antichains, which is often faster.
        """
        tau_list = self.tau_poset_elements
        num_tau = len(tau_list)
        
        # Precompute comparability matrix for the poset (τ, ⊆)
        comparable = [[False] * num_tau for _ in range(num_tau)]
        for i in range(num_tau):
            for j in range(i + 1, num_tau):
                if tau_list[i].issubset(tau_list[j]) or tau_list[j].issubset(tau_list[i]):
                    comparable[i][j] = comparable[j][i] = True

        memo = {}

        def count_antichains_recursive(k, current_antichain_indices_tuple):
            if k == num_tau:
                return 1
            
            state = (k, current_antichain_indices_tuple)
            if state in memo:
                return memo[state]

            # Case 1: element tau[k] is NOT in the antichain.
            count = count_antichains_recursive(k + 1, current_antichain_indices_tuple)
            
            # Case 2: element tau[k] IS in the antichain. Check for compatibility.
            is_valid_to_add = True
            current_indices = set(current_antichain_indices_tuple)
            for i in current_indices:
                if comparable[i][k]:
                    is_valid_to_add = False
                    break
            
            if is_valid_to_add:
                new_antichain_indices_tuple = current_antichain_indices_tuple + (k,)
                count += count_antichains_recursive(k + 1, new_antichain_indices_tuple)
            
            memo[state] = count
            return count

        return count_antichains_recursive(0, tuple())

# Assuming S is the set of divisors of 150
solver = PosetCounter(150)
result = solver.num_open_sets

print(f"Assuming S is the set of divisors of 150:")
print(f"The number of divisors is |S| = {len(solver.poset_elements)}")
print(f"The number of open sets in (D_S, τ) is |τ| = {len(solver.tau)}")
print(f"The number of open sets in P^-(D_S, τ) is the number of lower sets of (τ, ⊆).")
print(f"The final calculated number is: {result}")
print(f"The final equation is the result of this calculation. The number is {result}.")
