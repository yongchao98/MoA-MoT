import itertools
from collections import deque

def get_a_n_1d(n, search_span=10):
    """
    Calculates a(n) for the 1D case by exploring possible game evolutions.
    This is a brute-force approach that is only feasible for small n.
    """
    max_m_found = 1
    
    # We limit the search space for the initial configuration of n ones.
    # We place the first one at 0, then the rest within a certain span.
    initial_positions = range(1, search_span + n)
    
    for initial_config in itertools.combinations(initial_positions, n - 1):
        # A state is (m, occupied_cells). occupied_cells is a dict {pos: val}.
        # Start with the initial n ones.
        occupied = {0: 1}
        for pos in initial_config:
            occupied[pos] = 1

        queue = deque([(1, occupied)]) # (current_m, occupied_cells)
        visited_states = set()

        while queue:
            m, current_occupied = queue.popleft()
            
            # Use a frozenset of items for checking visited states
            frozen_state = frozenset(current_occupied.items())
            if frozen_state in visited_states:
                continue
            visited_states.add(frozen_state)

            if m > max_m_found:
                max_m_found = m
            
            # Find possible placements for k = m + 1
            k = m + 1
            candidates = []
            
            # Define search bounds for empty cells to check
            min_pos = min(current_occupied.keys())
            max_pos = max(current_occupied.keys())
            
            for p in range(min_pos - 2, max_pos + 3):
                if p not in current_occupied:
                    neighbor1_val = current_occupied.get(p - 1, 0)
                    neighbor2_val = current_occupied.get(p + 1, 0)
                    if neighbor1_val + neighbor2_val == k:
                        candidates.append(p)
            
            # If we found any valid placements, add them to the queue to explore
            if not candidates:
                # This branch of the game terminates
                if m > max_m_found:
                    max_m_found = m
            else:
                for cand_pos in candidates:
                    next_occupied = current_occupied.copy()
                    next_occupied[cand_pos] = k
                    queue.append((k, next_occupied))

    return max_m_found

def solve():
    """
    Provides the answers to the twelve questions.
    """
    # Answers to the first 9 (Yes/No) questions
    # 1. 3d, a(n)=inf? No
    # 2. 4d, a(n)=inf? No
    # 3. 5d, a(n)=inf? No
    # 4. 6d, a(n)=inf? No
    # 5. d>=1, a(n)<K*n? Yes
    # 6. 3d, a(n)>=9n-8? Yes
    # 7. 4d, a(n)>=17n-16? Yes
    # 8. 5d, a(n)<33n-32? Yes
    # 9. d>=2, a(n)>=(2^d+1)(n-1)+1? Yes
    answers_part1 = ["No", "No", "No", "No", "Yes", "Yes", "Yes", "Yes", "Yes"]
    
    # Answers to the last 3 (numerical) questions for the 1D case.
    # The code calculates a(2) and a(3). a(42) follows the same logic.
    a2 = get_a_n_1d(2)
    a3 = get_a_n_1d(3, search_span=6) # Smaller span for feasibility
    # Based on the analysis, the inability to form a 3 from 1s and 2s
    # is a structural problem in 1D that isn't resolved by adding more 1s.
    a42 = 2 
    answers_part2 = [a2, a3, a42]

    final_answers = answers_part1 + answers_part2
    print(",".join(map(str, final_answers)))

if __name__ == '__main__':
    solve()