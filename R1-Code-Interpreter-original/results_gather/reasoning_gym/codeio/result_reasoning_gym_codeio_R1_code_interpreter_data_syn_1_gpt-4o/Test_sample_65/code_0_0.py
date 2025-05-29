from typing import List
from collections import Counter

def main_solution(nums: List[int], quantity: List[int]) -> bool:
    counts = Counter(nums)
    left = sorted(counts.values())[-len(quantity):]
    quantity.sort(reverse=True)

    def can_fulfill_orders(left: List[int], quantity: List[int], customer: int) -> bool:
        if customer == len(quantity):
            return True
        for i in range(len(left)):
            if left[i] >= quantity[customer]:
                left[i] -= quantity[customer]
                if can_fulfill_orders(left, quantity, customer + 1):
                    return True
                left[i] += quantity[customer]
        return False

    return can_fulfill_orders(left, quantity, 0)

# Test case
nums = [1, 1, 1, 1]
quantity = [2, 2]
result = main_solution(nums, quantity)
print(result)